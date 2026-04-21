#!/usr/bin/env bash
CONFIG_YAML="$1"
PAIR_TSV="$2"
PROJECT_OUT="${3:-results}"
MAX_PARALLEL="${MAX_PARALLEL:-2}"

bash 01_submit_slurm_array.sh \
  --pipeline all \
  --pairs "$PAIR_TSV" \
  --results-dir "$PROJECT_OUT" \
  --max-parallel "$MAX_PARALLEL" \
  --config "$CONFIG_YAML"
