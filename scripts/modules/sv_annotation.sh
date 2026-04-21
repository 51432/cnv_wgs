#!/usr/bin/env bash
PAIR_TSV="$1"
MERGED_DIR="$2"
ASCAT_DIR="$3"
HPV_BREAKPOINTS="$4"
OUTDIR="$5"
TASK_ID="${SLURM_ARRAY_TASK_ID:-$6}"
SKIP_DONE="${SKIP_DONE:-true}"

source "$(dirname "$0")/../lib/common.sh"
ensure_log_dirs

row=$(get_pair_by_task_id "$PAIR_TSV" "$TASK_ID")
sample_id=$(echo "$row" | cut -f1)
sample_out="$OUTDIR/$sample_id"
ensure_dir "$sample_out"

if [ "$SKIP_DONE" = "true" ] && is_done "$sample_out/.done"; then
  echo "[sv_annotation] skip $sample_id"
  exit 0
fi

python3 "$(dirname "$0")/../python/sv_annotation.py" \
  --sample-id "$sample_id" \
  --merged-tsv "$MERGED_DIR/$sample_id/merged.sv.tsv" \
  --ascat-dir "$ASCAT_DIR/$sample_id" \
  --hpv-breakpoints "$HPV_BREAKPOINTS" \
  --outdir "$sample_out"
rc=$?
if [ $rc -ne 0 ]; then
  exit $rc
fi

mark_done "$sample_out/.done" "module=sv_annotation" "sample_id=$sample_id"
