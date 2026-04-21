#!/usr/bin/env bash

PAIRS="examples/pairbam.example.tsv"
GROUP=""
HPV="examples/hpv_breakpoints.example.tsv"
RESULTS_DIR="results"
MAX_PARALLEL=2
CONFIG="config/config.yaml.example"
PHASE="all"
START_STAGE=""
END_STAGE=""
SKIP_INPUT_CHECK=0

# repeatable: --max-parallel-stage stage=n
STAGE_PARALLEL_ARGS=()

usage() {
  cat <<'EOF'
Usage:
  bash 01_submit_slurm_array.sh \
    --pairs input/pairbam.tsv \
    --group input/group.tsv \
    --hpv-breakpoints input/hpv_breakpoints.tsv \
    --results-dir results_run \
    --phase precheck|sample|cohort|all \
    --start-stage ascat_prepare \
    --end-stage final_report \
    --max-parallel 2 \
    --max-parallel-stage sv_call_gridss=1 \
    --skip-input-check 0|1 \
    --config config/config.yaml

Phase:
  precheck: input_check
  sample:   soft_qc, ascat_prepare, ascat_run, sv_call_manta, sv_call_gridss,
            sv_postfilter, sv_merge, sv_annotation, hpv_link
  cohort:   cohort_summary, group_compare, final_report
  all:      full DAG (group.tsv 缺失时自动跳过 cohort 部分)
EOF
}

while [ $# -gt 0 ]; do
  case "$1" in
    --pairs) PAIRS="$2"; shift 2 ;;
    --group) GROUP="$2"; shift 2 ;;
    --hpv-breakpoints) HPV="$2"; shift 2 ;;
    --results-dir) RESULTS_DIR="$2"; shift 2 ;;
    --phase) PHASE="$2"; shift 2 ;;
    --max-parallel) MAX_PARALLEL="$2"; shift 2 ;;
    --max-parallel-stage) STAGE_PARALLEL_ARGS+=("$2"); shift 2 ;;
    --start-stage) START_STAGE="$2"; shift 2 ;;
    --end-stage) END_STAGE="$2"; shift 2 ;;
    --skip-input-check) SKIP_INPUT_CHECK="$2"; shift 2 ;;
    --config) CONFIG="$2"; shift 2 ;;
    -h|--help) usage; exit 0 ;;
    *) echo "Unknown arg: $1"; usage; exit 1 ;;
  esac
done

source scripts/lib/common.sh
ensure_log_dirs
ensure_results_layout "$RESULTS_DIR"

if [ ! -f "$PAIRS" ]; then
  echo "[error] pair file not found: $PAIRS"
  exit 1
fi
if [ ! -f "$CONFIG" ]; then
  echo "[error] config file not found: $CONFIG"
  exit 1
fi
snapshot_runtime_config "$CONFIG" "$PAIRS" "$RESULTS_DIR"

N=$(sample_count "$PAIRS")
if [ "$N" -lt 1 ]; then
  echo "[error] empty pair file"
  exit 1
fi

GROUP_EXISTS=0
GROUP_VALID_COLS=0
if [ -n "$GROUP" ] && [ -f "$GROUP" ]; then
  GROUP_EXISTS=1
  header=$(head -n 1 "$GROUP")
  ncol=$(echo "$header" | awk -F'\t' '{print NF}')
  if [ "$ncol" -gt 1 ]; then
    valid=$(echo "$header" | awk -F'\t' '{for(i=2;i<=NF;i++){if($i!=""){c++}}} END{print c+0}')
    if [ "$valid" -gt 0 ]; then
      GROUP_VALID_COLS=1
    fi
  fi
fi

all_stages=(input_check soft_qc ascat_prepare ascat_run sv_call_manta sv_call_gridss sv_postfilter sv_merge sv_annotation hpv_link cohort_summary group_compare final_report)
phase_precheck=(input_check)
phase_sample=(soft_qc ascat_prepare ascat_run sv_call_manta sv_call_gridss sv_postfilter sv_merge sv_annotation hpv_link)
phase_cohort=(cohort_summary group_compare final_report)

if [ "$PHASE" = "precheck" ]; then
  base=("${phase_precheck[@]}")
elif [ "$PHASE" = "sample" ]; then
  base=("${phase_sample[@]}")
elif [ "$PHASE" = "cohort" ]; then
  if [ "$GROUP_EXISTS" -ne 1 ]; then
    echo "[error] --phase cohort requires existing --group file"
    exit 1
  fi
  base=("${phase_cohort[@]}")
else
  base=("${all_stages[@]}")
fi

if [ -z "$START_STAGE" ]; then
  START_STAGE="${base[0]}"
fi
if [ -z "$END_STAGE" ]; then
  END_STAGE="${base[${#base[@]}-1]}"
fi

selected=()
active=0
for s in "${base[@]}"; do
  if [ "$s" = "$START_STAGE" ]; then active=1; fi
  if [ "$active" -eq 1 ]; then selected+=("$s"); fi
  if [ "$s" = "$END_STAGE" ]; then break; fi
done
if [ ${#selected[@]} -eq 0 ]; then
  echo "[error] invalid stage window for phase=$PHASE start=$START_STAGE end=$END_STAGE"
  exit 1
fi

# group 缺失时的 cohort 逻辑：all 阶段自动跳过 cohort 模块
if [ "$PHASE" = "all" ] && [ "$GROUP_EXISTS" -ne 1 ]; then
  echo "[warn] group.tsv missing: skip cohort_summary/group_compare/final_report, run precheck+sample only"
  filtered=()
  for s in "${selected[@]}"; do
    if [ "$s" = "cohort_summary" ] || [ "$s" = "group_compare" ] || [ "$s" = "final_report" ]; then
      continue
    fi
    filtered+=("$s")
  done
  selected=("${filtered[@]}")
fi

if [ ${#selected[@]} -eq 0 ]; then
  echo "[warn] no stages remain after optional-skip filtering, nothing to submit."
  exit 0
fi

declare -A STAGE_PARALLEL
for kv in "${STAGE_PARALLEL_ARGS[@]}"; do
  stage="${kv%%=*}"
  val="${kv##*=}"
  STAGE_PARALLEL[$stage]="$val"
done

stage_parallel() {
  st="$1"
  if [ -n "${STAGE_PARALLEL[$st]}" ]; then
    echo "${STAGE_PARALLEL[$st]}"
  else
    echo "$MAX_PARALLEL"
  fi
}

has_stage(){
  q="$1"
  for s in "${selected[@]}"; do
    if [ "$s" = "$q" ]; then return 0; fi
  done
  return 1
}

sample_done_exists() {
  dir="$1"
  missing=0
  while IFS=$'\t' read -r sid _; do
    if [ -z "$sid" ] || [ "$sid" = "sample_id" ]; then
      continue
    fi
    if [ ! -f "$dir/$sid/.done" ]; then
      echo "[error] missing upstream done: $dir/$sid/.done (required for rerun from middle stage)"
      missing=1
    fi
  done < "$PAIRS"
  [ "$missing" -eq 0 ]
}

require_stage_or_done() {
  stage="$1"
  upstream="$2"
  done_type="$3"
  done_path="$4"
  if has_stage "$stage" && ! has_stage "$upstream"; then
    if [ "$done_type" = "sample" ]; then
      sample_done_exists "$done_path" || exit 1
    else
      if [ ! -f "$done_path" ]; then
        echo "[error] missing upstream done: $done_path (required for rerun from middle stage)"
        exit 1
      fi
    fi
  fi
}

declare -A JOBID
submit_stage(){
  stage="$1"
  cmd="$2"
  deps="$3"
  dep_args=""
  if [ -n "$deps" ]; then
    dep_args="--dependency=afterok:$deps"
  fi
  jid=$(eval "sbatch $dep_args $cmd" | awk '{print $4}')
  JOBID[$stage]="$jid"
  echo "[submit] $stage => $jid deps=[$deps]"
}

dep_join(){
  out=""
  for d in "$@"; do
    if [ -n "$d" ]; then
      if [ -z "$out" ]; then out="$d"; else out="$out:$d"; fi
    fi
  done
  echo "$out"
}

echo "[info] phase=$PHASE selected stages=${selected[*]}"

# 中途重跑防护：若上游不在当前窗口，要求已有 done 标记
require_stage_or_done ascat_run ascat_prepare sample "$RESULTS_DIR/ascat_prepare"
require_stage_or_done sv_postfilter sv_call_manta sample "$RESULTS_DIR/sv/manta"
require_stage_or_done sv_postfilter sv_call_gridss sample "$RESULTS_DIR/sv/gridss"
require_stage_or_done sv_merge sv_postfilter sample "$RESULTS_DIR/sv/postfilter"
require_stage_or_done sv_annotation sv_merge sample "$RESULTS_DIR/sv/merged"
require_stage_or_done sv_annotation ascat_run sample "$RESULTS_DIR/ascat"
if has_stage hpv_link && [ -f "$HPV" ]; then
  require_stage_or_done hpv_link ascat_run sample "$RESULTS_DIR/ascat"
  require_stage_or_done hpv_link sv_annotation sample "$RESULTS_DIR/sv/annotation"
fi
require_stage_or_done cohort_summary ascat_run sample "$RESULTS_DIR/ascat"
require_stage_or_done cohort_summary sv_annotation sample "$RESULTS_DIR/sv/annotation"
require_stage_or_done group_compare cohort_summary single "$RESULTS_DIR/cohort/.done"
require_stage_or_done final_report cohort_summary single "$RESULTS_DIR/cohort/.done"

after_input=""
if [ "$SKIP_INPUT_CHECK" = "0" ]; then
  submit_stage input_check "slurm/input_check.sbatch $PAIRS $CONFIG $RESULTS_DIR/input_check" ""
  after_input="${JOBID[input_check]}"
else
  echo "[warn] skip input_check enabled; you must ensure upstream compatibility manually."
fi

if has_stage soft_qc; then
  p=$(stage_parallel soft_qc)
  submit_stage soft_qc "--array=1-$N%$p slurm/soft_qc.sbatch $PAIRS $CONFIG $RESULTS_DIR/soft_qc" "$after_input"
fi

if has_stage ascat_prepare; then
  p=$(stage_parallel ascat_prepare)
  submit_stage ascat_prepare "--array=1-$N%$p slurm/ascat_prepare.sbatch $PAIRS $CONFIG $RESULTS_DIR/ascat_prepare" "$after_input"
fi

if has_stage ascat_run; then
  p=$(stage_parallel ascat_run)
  d=$(dep_join "$after_input" "${JOBID[ascat_prepare]}")
  submit_stage ascat_run "--array=1-$N%$p slurm/ascat_run.sbatch $PAIRS $CONFIG $RESULTS_DIR/ascat_prepare $RESULTS_DIR/ascat" "$d"
fi

if has_stage sv_call_manta; then
  p=$(stage_parallel sv_call_manta)
  submit_stage sv_call_manta "--array=1-$N%$p slurm/sv_call_manta.sbatch $PAIRS $CONFIG $RESULTS_DIR/sv/manta" "$after_input"
fi

if has_stage sv_call_gridss; then
  p=$(stage_parallel sv_call_gridss)
  submit_stage sv_call_gridss "--array=1-$N%$p slurm/sv_call_gridss.sbatch $PAIRS $CONFIG $RESULTS_DIR/sv/gridss" "$after_input"
fi

if has_stage sv_postfilter; then
  p=$(stage_parallel sv_postfilter)
  d=$(dep_join "$after_input" "${JOBID[sv_call_manta]}" "${JOBID[sv_call_gridss]}")
  submit_stage sv_postfilter "--array=1-$N%$p slurm/sv_postfilter.sbatch $PAIRS $RESULTS_DIR/sv/manta $RESULTS_DIR/sv/gridss $RESULTS_DIR/sv/postfilter" "$d"
fi

if has_stage sv_merge; then
  p=$(stage_parallel sv_merge)
  d=$(dep_join "$after_input" "${JOBID[sv_postfilter]}")
  submit_stage sv_merge "--array=1-$N%$p slurm/sv_merge.sbatch $PAIRS $RESULTS_DIR/sv/postfilter $RESULTS_DIR/sv/merged" "$d"
fi

if has_stage sv_annotation; then
  p=$(stage_parallel sv_annotation)
  d=$(dep_join "$after_input" "${JOBID[sv_merge]}" "${JOBID[ascat_run]}")
  submit_stage sv_annotation "--array=1-$N%$p slurm/sv_annotation.sbatch $PAIRS $RESULTS_DIR/sv/merged $RESULTS_DIR/ascat $HPV $RESULTS_DIR/sv/annotation" "$d"
fi

if has_stage hpv_link; then
  if [ -f "$HPV" ]; then
    p=$(stage_parallel hpv_link)
    d=$(dep_join "$after_input" "${JOBID[ascat_run]}" "${JOBID[sv_annotation]}")
    submit_stage hpv_link "--array=1-$N%$p slurm/hpv_link.sbatch $PAIRS $HPV $RESULTS_DIR/ascat $RESULTS_DIR/sv/annotation $RESULTS_DIR/hpv_link" "$d"
  else
    echo "[auto-skip] hpv_link (missing hpv_breakpoints.tsv: $HPV)"
  fi
fi

if has_stage cohort_summary; then
  d=$(dep_join "$after_input" "${JOBID[ascat_run]}" "${JOBID[sv_annotation]}")
  submit_stage cohort_summary "slurm/cohort_summary.sbatch $RESULTS_DIR/ascat $RESULTS_DIR/sv/annotation $RESULTS_DIR/cohort" "$d"
fi

if has_stage group_compare; then
  if [ "$GROUP_EXISTS" -eq 1 ] && [ "$GROUP_VALID_COLS" -eq 1 ]; then
    d=$(dep_join "$after_input" "${JOBID[cohort_summary]}")
    if has_stage hpv_link; then
      d=$(dep_join "$d" "${JOBID[hpv_link]}")
    fi
    submit_stage group_compare "slurm/group_compare.sbatch $GROUP $RESULTS_DIR/cohort $RESULTS_DIR/group_compare" "$d"
  else
    echo "[auto-skip] group_compare (group.tsv missing or no valid group columns besides sample_id)"
  fi
fi

if has_stage final_report; then
  d=$(dep_join "$after_input" "${JOBID[cohort_summary]}" "${JOBID[soft_qc]}" "${JOBID[group_compare]}" "${JOBID[hpv_link]}")
  submit_stage final_report "slurm/final_report.sbatch $RESULTS_DIR $RESULTS_DIR/final_report" "$d"
fi
