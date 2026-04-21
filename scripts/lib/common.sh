#!/usr/bin/env bash

PROJECT_ROOT=$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)

ensure_dir() {
  d="$1"
  if [ ! -d "$d" ]; then
    mkdir -p "$d"
    echo "[mkdir] $d"
  else
    echo "[keep]  $d"
  fi
}

ensure_log_dirs() {
  ensure_dir "$PROJECT_ROOT/logs"
  for d in input_check soft_qc ascat_prepare ascat manta gridss sv_postfilter sv_merge sv_annotation cohort_summary group_compare hpv_link report snv_indel_hook submit; do
    ensure_dir "$PROJECT_ROOT/logs/$d"
  done
}

ensure_results_layout() {
  root="$1"
  ensure_dir "$root"
  ensure_dir "$root/00.config"
  ensure_dir "$root/input_check"
  ensure_dir "$root/soft_qc"
  ensure_dir "$root/ascat_prepare"
  ensure_dir "$root/ascat"
  ensure_dir "$root/sv"
  ensure_dir "$root/sv/manta"
  ensure_dir "$root/sv/gridss"
  ensure_dir "$root/sv/postfilter"
  ensure_dir "$root/sv/merged"
  ensure_dir "$root/sv/annotation"
  ensure_dir "$root/cohort"
  ensure_dir "$root/group_compare"
  ensure_dir "$root/hpv_link"
  ensure_dir "$root/final_report"
}

snapshot_runtime_config() {
  cfg="$1"
  pairs="$2"
  outdir="$3/00.config"
  ensure_dir "$outdir"
  ts=$(date '+%Y%m%d_%H%M%S')
  cp "$cfg" "$outdir/config.${ts}.yaml"
  cp "$pairs" "$outdir/pairbam.${ts}.tsv"
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
  stage_val=""
  sample_val=""
  echo "status=done" > "$done_file"
  echo "timestamp=$ts" >> "$done_file"
  for kv in "$@"; do
    echo "$kv" >> "$done_file"
    case "$kv" in
      module=*) stage_val="${kv#module=}" ;;
      stage=*) stage_val="${kv#stage=}" ;;
      sample_id=*) sample_val="${kv#sample_id=}" ;;
    esac
  done
  if [ -n "$stage_val" ]; then
    echo "stage=$stage_val" >> "$done_file"
  fi
  if [ -n "$sample_val" ]; then
    echo "sample_id=$sample_val" >> "$done_file"
  fi
}
