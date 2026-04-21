#!/usr/bin/env bash
PROJECT_OUTDIR="$1"
OUTDIR="$2"
source "$(dirname "$0")/../lib/common.sh"
ensure_log_dirs
ensure_dir "$OUTDIR"
python3 "$(dirname "$0")/../python/final_report.py" --project-outdir "$PROJECT_OUTDIR" --outdir "$OUTDIR"
rc=$?
if [ $rc -ne 0 ]; then exit $rc; fi
mark_done "$OUTDIR/.done" "module=final_report"
