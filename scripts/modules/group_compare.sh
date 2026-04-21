#!/usr/bin/env bash
GROUP_TSV="$1"
COHORT_DIR="$2"
OUTDIR="$3"
source "$(dirname "$0")/../lib/common.sh"
ensure_log_dirs
ensure_dir "$OUTDIR"
python3 "$(dirname "$0")/../python/group_compare.py" --group-tsv "$GROUP_TSV" --cohort-dir "$COHORT_DIR" --outdir "$OUTDIR"
rc=$?
if [ $rc -ne 0 ]; then exit $rc; fi
mark_done "$OUTDIR/.done" "module=group_compare"
