#!/usr/bin/env bash
CONFIG_YAML="$1"
PAIR_TSV="$2"
PROJECT_OUT="${3:-results}"

source "$(dirname "$0")/lib/common.sh"
ensure_log_dirs
ensure_dir "$PROJECT_OUT"

# 1) input check
sbatch slurm/input_check.sbatch "$PAIR_TSV" "$CONFIG_YAML" "$PROJECT_OUT/input_check"

# 2) 样本级模块建议按依赖逐步提交（可用 --dependency=afterok）
N=$(sample_count "$PAIR_TSV")
P="${MAX_PARALLEL:-2}"
sbatch --array=1-"$N"%"$P" slurm/soft_qc.sbatch "$PAIR_TSV" "$PROJECT_OUT/soft_qc"
sbatch --array=1-"$N"%"$P" slurm/ascat_prepare.sbatch "$PAIR_TSV" "$CONFIG_YAML" "$PROJECT_OUT/ascat_prepare"
sbatch --array=1-"$N"%"$P" slurm/ascat_run.sbatch "$PAIR_TSV" "$CONFIG_YAML" "$PROJECT_OUT/ascat_prepare" "$PROJECT_OUT/ascat"
sbatch --array=1-"$N"%"$P" slurm/sv_call_manta.sbatch "$PAIR_TSV" "$CONFIG_YAML" "$PROJECT_OUT/sv/manta"
sbatch --array=1-"$N"%"$P" slurm/sv_call_gridss.sbatch "$PAIR_TSV" "$CONFIG_YAML" "$PROJECT_OUT/sv/gridss"
sbatch --array=1-"$N"%"$P" slurm/sv_postfilter.sbatch "$PAIR_TSV" "$PROJECT_OUT/sv/manta" "$PROJECT_OUT/sv/gridss" "$PROJECT_OUT/sv/postfilter"

# 3) 队列级
sbatch slurm/sv_merge.sbatch "$PAIR_TSV" "$PROJECT_OUT/sv/postfilter" "$PROJECT_OUT/sv/merged"
sbatch slurm/sv_annotation.sbatch "$PROJECT_OUT/sv/merged" "$PROJECT_OUT/ascat" "${HPV_BREAKPOINTS:-examples/hpv_breakpoints.example.tsv}" "$PROJECT_OUT/sv/annotation"
sbatch slurm/cohort_summary.sbatch "$PROJECT_OUT/ascat" "$PROJECT_OUT/sv/annotation" "$PROJECT_OUT/cohort"
