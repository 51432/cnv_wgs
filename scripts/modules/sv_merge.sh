#!/usr/bin/env bash
PAIR_TSV="$1"
POSTFILTER_DIR="$2"
OUTDIR="$3"
source "$(dirname "$0")/../lib/common.sh"
ensure_log_dirs
ensure_dir "$OUTDIR"
python3 "$(dirname "$0")/../python/sv_merge.py" --pairbam "$PAIR_TSV" --postfilter-dir "$POSTFILTER_DIR" --outdir "$OUTDIR"
rc=$?
if [ $rc -ne 0 ]; then exit $rc; fi
mark_done "$OUTDIR/.done" "module=sv_merge"
