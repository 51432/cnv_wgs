#!/usr/bin/env bash
ASCAT_DIR="$1"
SV_ANNOT_DIR="$2"
OUTDIR="$3"
source "$(dirname "$0")/../lib/common.sh"
ensure_log_dirs
ensure_dir "$OUTDIR"
python3 "$(dirname "$0")/../python/cohort_summary.py" --ascat-dir "$ASCAT_DIR" --sv-annot-dir "$SV_ANNOT_DIR" --outdir "$OUTDIR"
rc=$?
if [ $rc -ne 0 ]; then exit $rc; fi
mark_done "$OUTDIR/.done" "module=cohort_summary"
