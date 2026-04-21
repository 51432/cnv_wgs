#!/usr/bin/env bash
HPV_BREAKPOINTS="$1"
ASCAT_DIR="$2"
SV_ANNOT_DIR="$3"
OUTDIR="$4"
source "$(dirname "$0")/../lib/common.sh"
ensure_log_dirs
ensure_dir "$OUTDIR"
python3 "$(dirname "$0")/../python/hpv_link.py" --hpv-breakpoints "$HPV_BREAKPOINTS" --ascat-dir "$ASCAT_DIR" --sv-annot-dir "$SV_ANNOT_DIR" --outdir "$OUTDIR"
rc=$?
if [ $rc -ne 0 ]; then exit $rc; fi
mark_done "$OUTDIR/.done" "module=hpv_link"
