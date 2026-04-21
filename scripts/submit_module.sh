#!/usr/bin/env bash

MODULE="$1"
CONFIG_YAML="$2"
PAIR_TSV="$3"
EXTRA_ARGS="$4"

source "$(dirname "$0")/lib/common.sh"
ensure_log_dirs

case "$MODULE" in
  soft_qc|ascat_prepare|ascat_run|sv_call_manta|sv_call_gridss|sv_postfilter)
    N=$(sample_count "$PAIR_TSV")
    P="${MAX_PARALLEL:-2}"
    sbatch --array=1-"$N"%"$P" "$EXTRA_ARGS"
    ;;
  *)
    sbatch "$EXTRA_ARGS"
    ;;
esac
