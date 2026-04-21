#!/usr/bin/env bash
PAIR_TSV="$1"
CONFIG_YAML="$2"
OUTDIR="$3"
TASK_ID="${SLURM_ARRAY_TASK_ID:-$4}"
SKIP_DONE="${SKIP_DONE:-true}"
source "$(dirname "$0")/../lib/common.sh"
ensure_log_dirs
row=$(get_pair_by_task_id "$PAIR_TSV" "$TASK_ID")
sample_id=$(echo "$row" | cut -f1)
tumor_bam=$(echo "$row" | cut -f2)
normal_bam=$(echo "$row" | cut -f3)
sample_dir="$OUTDIR/$sample_id"
ensure_dir "$sample_dir"
if [ "$SKIP_DONE" = "true" ] && is_done "$sample_dir/.done"; then echo "[sv_call_manta] skip $sample_id"; exit 0; fi
python3 "$(dirname "$0")/../python/sv_call_manta.py" --sample-id "$sample_id" --tumor-bam "$tumor_bam" --normal-bam "$normal_bam" --config "$CONFIG_YAML" --outdir "$sample_dir"
rc=$?
if [ $rc -ne 0 ]; then exit $rc; fi
mark_done "$sample_dir/.done" "module=sv_call_manta" "sample_id=$sample_id"
