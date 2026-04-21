#!/usr/bin/env bash
MERGED_DIR="$1"
ASCAT_DIR="$2"
HPV_BREAKPOINTS="$3"
OUTDIR="$4"
source "$(dirname "$0")/../lib/common.sh"
ensure_log_dirs
ensure_dir "$OUTDIR"
python3 "$(dirname "$0")/../python/sv_annotation.py" --merged-dir "$MERGED_DIR" --ascat-dir "$ASCAT_DIR" --hpv-breakpoints "$HPV_BREAKPOINTS" --outdir "$OUTDIR"
rc=$?
if [ $rc -ne 0 ]; then exit $rc; fi
mark_done "$OUTDIR/.done" "module=sv_annotation"
