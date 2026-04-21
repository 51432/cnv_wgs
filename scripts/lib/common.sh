#!/usr/bin/env bash

PROJECT_ROOT=$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)

ensure_dir() {
  d="$1"
  if [ ! -d "$d" ]; then
    mkdir -p "$d"
  fi
}

ensure_log_dirs() {
  ensure_dir "$PROJECT_ROOT/logs"
  for d in input_check soft_qc ascat_prepare ascat manta gridss sv_postfilter sv_merge sv_annotation cohort_summary group_compare hpv_link report snv_indel_hook; do
    ensure_dir "$PROJECT_ROOT/logs/$d"
  done
}

sample_count() {
  pair_tsv="$1"
  awk 'NR>1{n++} END{print n+0}' "$pair_tsv"
}

get_pair_by_task_id() {
  pair_tsv="$1"
  task_id="$2"
  awk -F'\t' -v id="$task_id" 'NR==1{next} NR==(id+1){print $1"\t"$2"\t"$3}' "$pair_tsv"
}

is_done() {
  done_file="$1"
  if [ -f "$done_file" ]; then
    return 0
  fi
  return 1
}

mark_done() {
  done_file="$1"
  shift
  ts=$(date '+%Y-%m-%d %H:%M:%S')
  echo "status=done" > "$done_file"
  echo "timestamp=$ts" >> "$done_file"
  for kv in "$@"; do
    echo "$kv" >> "$done_file"
  done
}

