#!/usr/bin/env bash

PIPELINE="all"
PAIRS="examples/pairbam.example.tsv"
GROUP="examples/group.example.tsv"
HPV="examples/hpv_breakpoints.example.tsv"
MODE="wgs"
RESULTS_DIR="results"
MAX_PARALLEL=2
START_STAGE="input_check"
END_STAGE="final_report"
ENABLE_CONTAMINATION=0
ENABLE_ORIENTATION=0
ENABLE_ANNOTATION=1
CONFIG="config/config.yaml.example"

usage() {
  cat <<'EOF'
Usage:
  bash 01_submit_slurm_array.sh \
    --pipeline all|phase1|phase2|phase3 \
    --pairs input/pairbam.tsv \
    --group input/group.tsv \
    --hpv-breakpoints input/hpv_breakpoints.tsv \
    --mode wgs|wes \
    --results-dir results \
    --max-parallel 2 \
    --start-stage input_check \
    --end-stage final_report \
    --enable-contamination 0|1 \
    --enable-orientation 0|1 \
    --enable-annotation 0|1 \
    --config config/config.yaml
EOF
}

while [ $# -gt 0 ]; do
  case "$1" in
    --pipeline) PIPELINE="$2"; shift 2 ;;
    --pairs) PAIRS="$2"; shift 2 ;;
    --group) GROUP="$2"; shift 2 ;;
    --hpv-breakpoints) HPV="$2"; shift 2 ;;
    --mode) MODE="$2"; shift 2 ;;
    --results-dir) RESULTS_DIR="$2"; shift 2 ;;
    --max-parallel) MAX_PARALLEL="$2"; shift 2 ;;
    --start-stage) START_STAGE="$2"; shift 2 ;;
    --end-stage) END_STAGE="$2"; shift 2 ;;
    --enable-contamination) ENABLE_CONTAMINATION="$2"; shift 2 ;;
    --enable-orientation) ENABLE_ORIENTATION="$2"; shift 2 ;;
    --enable-annotation) ENABLE_ANNOTATION="$2"; shift 2 ;;
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
N=$(sample_count "$PAIRS")
if [ "$N" -lt 1 ]; then
  echo "[error] empty pair file"
  exit 1
fi

all_stages=(input_check soft_qc ascat_prepare ascat_run sv_call_manta sv_call_gridss sv_postfilter sv_merge sv_annotation hpv_link cohort_summary group_compare final_report)
phase1=(input_check soft_qc ascat_prepare ascat_run sv_call_manta sv_call_gridss sv_postfilter sv_merge sv_annotation hpv_link)
phase2=(cohort_summary group_compare)
phase3=(final_report)

normalize_stage(){
  if [ "$1" = "filter" ]; then echo "sv_postfilter"; else echo "$1"; fi
}
START_STAGE=$(normalize_stage "$START_STAGE")
END_STAGE=$(normalize_stage "$END_STAGE")

if [ "$PIPELINE" = "phase1" ]; then
  base=("${phase1[@]}")
elif [ "$PIPELINE" = "phase2" ]; then
  base=("${phase2[@]}")
elif [ "$PIPELINE" = "phase3" ]; then
  base=("${phase3[@]}")
else
  base=("${all_stages[@]}")
fi

selected=()
active=0
for s in "${base[@]}"; do
  if [ "$s" = "$START_STAGE" ]; then active=1; fi
  if [ "$active" -eq 1 ]; then selected+=("$s"); fi
  if [ "$s" = "$END_STAGE" ]; then break; fi
done

if [ ${#selected[@]} -eq 0 ]; then
  echo "[error] invalid stage window"
  exit 1
fi

has_stage(){
  q="$1"
  for s in "${selected[@]}"; do
    if [ "$s" = "$q" ]; then return 0; fi
  done
  return 1
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

echo "[info] mode=$MODE max_parallel=$MAX_PARALLEL stages=${selected[*]}"

after_input=""
if has_stage input_check; then
  submit_stage input_check "slurm/input_check.sbatch $PAIRS $CONFIG $RESULTS_DIR/input_check" ""
  after_input="${JOBID[input_check]}"
fi

if has_stage soft_qc; then
  submit_stage soft_qc "--array=1-$N%$MAX_PARALLEL slurm/soft_qc.sbatch $PAIRS $RESULTS_DIR/soft_qc" "$after_input"
fi

if has_stage ascat_prepare; then
  submit_stage ascat_prepare "--array=1-$N%$MAX_PARALLEL slurm/ascat_prepare.sbatch $PAIRS $CONFIG $RESULTS_DIR/ascat_prepare" "$after_input"
fi

if has_stage ascat_run; then
  d="${JOBID[ascat_prepare]}"
  submit_stage ascat_run "--array=1-$N%$MAX_PARALLEL slurm/ascat_run.sbatch $PAIRS $CONFIG $RESULTS_DIR/ascat_prepare $RESULTS_DIR/ascat" "$d"
fi

if has_stage sv_call_manta; then
  submit_stage sv_call_manta "--array=1-$N%$MAX_PARALLEL slurm/sv_call_manta.sbatch $PAIRS $CONFIG $RESULTS_DIR/sv/manta" "$after_input"
fi

if has_stage sv_call_gridss; then
  submit_stage sv_call_gridss "--array=1-$N%$MAX_PARALLEL slurm/sv_call_gridss.sbatch $PAIRS $CONFIG $RESULTS_DIR/sv/gridss" "$after_input"
fi

if has_stage sv_postfilter; then
  d=$(dep_join "${JOBID[sv_call_manta]}" "${JOBID[sv_call_gridss]}")
  submit_stage sv_postfilter "--array=1-$N%$MAX_PARALLEL slurm/sv_postfilter.sbatch $PAIRS $RESULTS_DIR/sv/manta $RESULTS_DIR/sv/gridss $RESULTS_DIR/sv/postfilter" "$d"
fi

if has_stage sv_merge; then
  d="${JOBID[sv_postfilter]}"
  submit_stage sv_merge "--array=1-$N%$MAX_PARALLEL slurm/sv_merge.sbatch $PAIRS $RESULTS_DIR/sv/postfilter $RESULTS_DIR/sv/merged" "$d"
fi

if has_stage sv_annotation; then
  if [ "$ENABLE_ANNOTATION" = "1" ]; then
    d=$(dep_join "${JOBID[sv_merge]}" "${JOBID[ascat_run]}")
    submit_stage sv_annotation "--array=1-$N%$MAX_PARALLEL slurm/sv_annotation.sbatch $PAIRS $RESULTS_DIR/sv/merged $RESULTS_DIR/ascat $HPV $RESULTS_DIR/sv/annotation" "$d"
  else
    echo "[skip] sv_annotation (enable-annotation=0)"
  fi
fi

if has_stage hpv_link; then
  if [ -f "$HPV" ]; then
    d=$(dep_join "${JOBID[ascat_run]}" "${JOBID[sv_annotation]}")
    submit_stage hpv_link "--array=1-$N%$MAX_PARALLEL slurm/hpv_link.sbatch $PAIRS $HPV $RESULTS_DIR/ascat $RESULTS_DIR/sv/annotation $RESULTS_DIR/hpv_link" "$d"
  else
    echo "[skip] hpv_link (missing hpv breakpoints)"
  fi
fi

if has_stage cohort_summary; then
  d=$(dep_join "${JOBID[ascat_run]}" "${JOBID[sv_annotation]}" "${JOBID[hpv_link]}")
  submit_stage cohort_summary "slurm/cohort_summary.sbatch $RESULTS_DIR/ascat $RESULTS_DIR/sv/annotation $RESULTS_DIR/cohort" "$d"
fi

if has_stage group_compare; then
  if [ -f "$GROUP" ]; then
    d="${JOBID[cohort_summary]}"
    if has_stage hpv_link; then
      d=$(dep_join "$d" "${JOBID[hpv_link]}")
    fi
    submit_stage group_compare "slurm/group_compare.sbatch $GROUP $RESULTS_DIR/cohort $RESULTS_DIR/group_compare" "$d"
  else
    echo "[skip] group_compare (missing group.tsv)"
  fi
fi

if has_stage final_report; then
  d=$(dep_join "${JOBID[cohort_summary]}" "${JOBID[group_compare]}" "${JOBID[hpv_link]}")
  submit_stage final_report "slurm/final_report.sbatch $RESULTS_DIR $RESULTS_DIR/final_report" "$d"
fi

