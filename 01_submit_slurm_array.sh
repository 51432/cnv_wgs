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

Stage list:
  input_check, soft_qc, ascat_prepare, ascat_run,
  sv_call_manta, sv_call_gridss, sv_postfilter,
  sv_merge, sv_annotation, cohort_summary,
  group_compare, hpv_link, final_report
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
  echo "[error] pair file has no sample rows: $PAIRS"
  exit 1
fi

all_stages=(input_check soft_qc ascat_prepare ascat_run sv_call_manta sv_call_gridss sv_postfilter sv_merge sv_annotation cohort_summary group_compare hpv_link final_report)

phase1=(input_check soft_qc ascat_prepare ascat_run sv_call_manta sv_call_gridss sv_postfilter)
phase2=(sv_merge sv_annotation cohort_summary)
phase3=(group_compare hpv_link final_report)

if [ "$PIPELINE" = "phase1" ]; then
  stages=("${phase1[@]}")
elif [ "$PIPELINE" = "phase2" ]; then
  stages=("${phase2[@]}")
elif [ "$PIPELINE" = "phase3" ]; then
  stages=("${phase3[@]}")
else
  stages=("${all_stages[@]}")
fi

normalize_stage() {
  s="$1"
  if [ "$s" = "filter" ]; then
    echo "sv_postfilter"
    return
  fi
  echo "$s"
}

START_STAGE=$(normalize_stage "$START_STAGE")
END_STAGE=$(normalize_stage "$END_STAGE")

subset=()
active=0
for st in "${stages[@]}"; do
  if [ "$st" = "$START_STAGE" ]; then
    active=1
  fi
  if [ "$active" -eq 1 ]; then
    subset+=("$st")
  fi
  if [ "$st" = "$END_STAGE" ]; then
    break
  fi
done

if [ ${#subset[@]} -eq 0 ]; then
  echo "[error] stage window invalid: start=$START_STAGE end=$END_STAGE pipeline=$PIPELINE"
  exit 1
fi

echo "[info] mode=$MODE pipeline=$PIPELINE stages=${subset[*]} max_parallel=$MAX_PARALLEL"
echo "[info] contamination=$ENABLE_CONTAMINATION orientation=$ENABLE_ORIENTATION annotation=$ENABLE_ANNOTATION"

prev_job=""
submit_cmd() {
  stage="$1"
  cmd="$2"
  if [ -n "$prev_job" ]; then
    jid=$(eval "$cmd --dependency=afterok:$prev_job" | awk '{print $4}')
  else
    jid=$(eval "$cmd" | awk '{print $4}')
  fi
  echo "[submit] $stage => $jid"
  prev_job="$jid"
}

for stage in "${subset[@]}"; do
  case "$stage" in
    input_check)
      submit_cmd "$stage" "sbatch slurm/input_check.sbatch $PAIRS $CONFIG $RESULTS_DIR/input_check"
      ;;
    soft_qc)
      submit_cmd "$stage" "sbatch --array=1-$N%$MAX_PARALLEL slurm/soft_qc.sbatch $PAIRS $RESULTS_DIR/soft_qc"
      ;;
    ascat_prepare)
      submit_cmd "$stage" "sbatch --array=1-$N%$MAX_PARALLEL slurm/ascat_prepare.sbatch $PAIRS $CONFIG $RESULTS_DIR/ascat_prepare"
      ;;
    ascat_run)
      submit_cmd "$stage" "sbatch --array=1-$N%$MAX_PARALLEL slurm/ascat_run.sbatch $PAIRS $CONFIG $RESULTS_DIR/ascat_prepare $RESULTS_DIR/ascat"
      ;;
    sv_call_manta)
      submit_cmd "$stage" "sbatch --array=1-$N%$MAX_PARALLEL slurm/sv_call_manta.sbatch $PAIRS $CONFIG $RESULTS_DIR/sv/manta"
      ;;
    sv_call_gridss)
      submit_cmd "$stage" "sbatch --array=1-$N%$MAX_PARALLEL slurm/sv_call_gridss.sbatch $PAIRS $CONFIG $RESULTS_DIR/sv/gridss"
      ;;
    sv_postfilter)
      submit_cmd "$stage" "sbatch --array=1-$N%$MAX_PARALLEL slurm/sv_postfilter.sbatch $PAIRS $RESULTS_DIR/sv/manta $RESULTS_DIR/sv/gridss $RESULTS_DIR/sv/postfilter"
      ;;
    sv_merge)
      submit_cmd "$stage" "sbatch slurm/sv_merge.sbatch $PAIRS $RESULTS_DIR/sv/postfilter $RESULTS_DIR/sv/merged"
      ;;
    sv_annotation)
      if [ "$ENABLE_ANNOTATION" = "1" ]; then
        submit_cmd "$stage" "sbatch slurm/sv_annotation.sbatch $RESULTS_DIR/sv/merged $RESULTS_DIR/ascat $HPV $RESULTS_DIR/sv/annotation"
      else
        echo "[skip] $stage (enable-annotation=0)"
      fi
      ;;
    cohort_summary)
      submit_cmd "$stage" "sbatch slurm/cohort_summary.sbatch $RESULTS_DIR/ascat $RESULTS_DIR/sv/annotation $RESULTS_DIR/cohort"
      ;;
    group_compare)
      if [ -f "$GROUP" ]; then
        submit_cmd "$stage" "sbatch slurm/group_compare.sbatch $GROUP $RESULTS_DIR/cohort $RESULTS_DIR/group_compare"
      else
        echo "[skip] $stage (group file not found: $GROUP)"
      fi
      ;;
    hpv_link)
      if [ -f "$HPV" ]; then
        submit_cmd "$stage" "sbatch slurm/hpv_link.sbatch $HPV $RESULTS_DIR/ascat $RESULTS_DIR/sv/annotation $RESULTS_DIR/hpv_link"
      else
        echo "[skip] $stage (hpv breakpoints not found: $HPV)"
      fi
      ;;
    final_report)
      submit_cmd "$stage" "sbatch slurm/final_report.sbatch $RESULTS_DIR $RESULTS_DIR/final_report"
      ;;
    *)
      echo "[warn] unknown stage: $stage"
      ;;
  esac
done
