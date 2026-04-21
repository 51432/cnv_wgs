#!/usr/bin/env bash

PAIR_TSV="$1"
CONFIG_YAML="$2"
OUTDIR="$3"

source "$(dirname "$0")/../lib/common.sh"
ensure_log_dirs
ensure_dir "$OUTDIR"

python3 "$(dirname "$0")/../python/input_check.py" \
  --pairbam "$PAIR_TSV" \
  --config "$CONFIG_YAML" \
  --outdir "$OUTDIR"
rc=$?
if [ $rc -ne 0 ]; then
  echo "[input_check] failed"
  exit $rc
fi

mark_done "$OUTDIR/.done" "module=input_check" "pairbam=$PAIR_TSV"
