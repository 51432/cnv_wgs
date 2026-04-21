#!/usr/bin/env bash
PAIR_TSV="$1"
MANTA_DIR="$2"
GRIDSS_DIR="$3"
OUTDIR="$4"
TASK_ID="${SLURM_ARRAY_TASK_ID:-$5}"
SKIP_DONE="${SKIP_DONE:-true}"
source "$(dirname "$0")/../lib/common.sh"
ensure_log_dirs
row=$(get_pair_by_task_id "$PAIR_TSV" "$TASK_ID")
sample_id=$(echo "$row" | cut -f1)
sample_out="$OUTDIR/$sample_id"
ensure_dir "$sample_out"
if [ "$SKIP_DONE" = "true" ] && is_done "$sample_out/.done"; then echo "[sv_postfilter] skip $sample_id"; exit 0; fi
python3 "$(dirname "$0")/../python/sv_postfilter.py" --sample-id "$sample_id" --manta-vcf "$MANTA_DIR/$sample_id/somaticSV.vcf.gz" --gridss-vcf "$GRIDSS_DIR/$sample_id/gridss.somatic.vcf.gz" --outdir "$sample_out"
rc=$?
if [ $rc -ne 0 ]; then exit $rc; fi
mark_done "$sample_out/.done" "module=sv_postfilter" "sample_id=$sample_id"
