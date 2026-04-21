#!/usr/bin/env bash
PAIR_TSV="$1"
OUTDIR="$2"
source "$(dirname "$0")/../lib/common.sh"
ensure_log_dirs
ensure_dir "$OUTDIR"
echo "snv_indel_hook is reserved for future SNV/indel integration" > "$OUTDIR/README.txt"
mark_done "$OUTDIR/.done" "module=snv_indel_hook"
