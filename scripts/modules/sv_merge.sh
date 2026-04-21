#!/usr/bin/env bash
PAIR_TSV="$1"
POSTFILTER_DIR="$2"
OUTDIR="$3"
TASK_ID="${SLURM_ARRAY_TASK_ID:-$4}"
SKIP_DONE="${SKIP_DONE:-true}"

source "$(dirname "$0")/../lib/common.sh"
ensure_log_dirs

row=$(get_pair_by_task_id "$PAIR_TSV" "$TASK_ID")
sample_id=$(echo "$row" | cut -f1)
sample_out="$OUTDIR/$sample_id"
ensure_dir "$sample_out"

if [ "$SKIP_DONE" = "true" ] && is_done "$sample_out/.done"; then
  echo "[sv_merge] skip $sample_id"
  exit 0
fi

python3 "$(dirname "$0")/../python/sv_merge.py" \
  --sample-id "$sample_id" \
  --postfilter-dir "$POSTFILTER_DIR/$sample_id" \
  --outdir "$sample_out"
rc=$?
if [ $rc -ne 0 ]; then
  exit $rc
fi

mark_done "$sample_out/.done" "module=sv_merge" "sample_id=$sample_id"
