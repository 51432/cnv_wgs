#!/usr/bin/env bash
PAIR_TSV="$1"
CONFIG_YAML="$2"
PREP_DIR="$3"
OUTDIR="$4"
TASK_ID="${SLURM_ARRAY_TASK_ID:-$5}"
SKIP_DONE="${SKIP_DONE:-true}"
source "$(dirname "$0")/../lib/common.sh"
ensure_log_dirs
row=$(get_pair_by_task_id "$PAIR_TSV" "$TASK_ID")
sample_id=$(echo "$row" | cut -f1)
sample_prep="$PREP_DIR/$sample_id"
sample_out="$OUTDIR/$sample_id"
ensure_dir "$sample_out"
if [ "$SKIP_DONE" = "true" ] && is_done "$sample_out/.done"; then echo "[ascat_run] skip $sample_id"; exit 0; fi
Rscript "$(dirname "$0")/../r/ascat_run.R" --sample_id "$sample_id" --prep_dir "$sample_prep" --config "$CONFIG_YAML" --outdir "$sample_out"
rc=$?
if [ $rc -ne 0 ]; then exit $rc; fi
mark_done "$sample_out/.done" "module=ascat_run" "sample_id=$sample_id"
