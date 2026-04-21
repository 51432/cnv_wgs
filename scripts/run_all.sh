#!/usr/bin/env bash
CONFIG_YAML="$1"
PAIR_TSV="$2"
PROJECT_OUT="${3:-results}"

source "$(dirname "$0")/lib/common.sh"
ensure_log_dirs
ensure_results_layout "$PROJECT_OUT"

N=$(sample_count "$PAIR_TSV")
P="${MAX_PARALLEL:-2}"

jid1=$(sbatch slurm/input_check.sbatch "$PAIR_TSV" "$CONFIG_YAML" "$PROJECT_OUT/input_check" | awk '{print $4}')
jid2=$(sbatch --dependency=afterok:"$jid1" --array=1-"$N"%"$P" slurm/soft_qc.sbatch "$PAIR_TSV" "$PROJECT_OUT/soft_qc" | awk '{print $4}')
jid3=$(sbatch --dependency=afterok:"$jid2" --array=1-"$N"%"$P" slurm/ascat_prepare.sbatch "$PAIR_TSV" "$CONFIG_YAML" "$PROJECT_OUT/ascat_prepare" | awk '{print $4}')
jid4=$(sbatch --dependency=afterok:"$jid3" --array=1-"$N"%"$P" slurm/ascat_run.sbatch "$PAIR_TSV" "$CONFIG_YAML" "$PROJECT_OUT/ascat_prepare" "$PROJECT_OUT/ascat" | awk '{print $4}')
jid5=$(sbatch --dependency=afterok:"$jid4" --array=1-"$N"%"$P" slurm/sv_call_manta.sbatch "$PAIR_TSV" "$CONFIG_YAML" "$PROJECT_OUT/sv/manta" | awk '{print $4}')
jid6=$(sbatch --dependency=afterok:"$jid5" --array=1-"$N"%"$P" slurm/sv_call_gridss.sbatch "$PAIR_TSV" "$CONFIG_YAML" "$PROJECT_OUT/sv/gridss" | awk '{print $4}')
jid7=$(sbatch --dependency=afterok:"$jid6" --array=1-"$N"%"$P" slurm/sv_postfilter.sbatch "$PAIR_TSV" "$PROJECT_OUT/sv/manta" "$PROJECT_OUT/sv/gridss" "$PROJECT_OUT/sv/postfilter" | awk '{print $4}')
jid8=$(sbatch --dependency=afterok:"$jid7" slurm/sv_merge.sbatch "$PAIR_TSV" "$PROJECT_OUT/sv/postfilter" "$PROJECT_OUT/sv/merged" | awk '{print $4}')
jid9=$(sbatch --dependency=afterok:"$jid8" slurm/sv_annotation.sbatch "$PROJECT_OUT/sv/merged" "$PROJECT_OUT/ascat" "${HPV_BREAKPOINTS:-examples/hpv_breakpoints.example.tsv}" "$PROJECT_OUT/sv/annotation" | awk '{print $4}')
jid10=$(sbatch --dependency=afterok:"$jid9" slurm/cohort_summary.sbatch "$PROJECT_OUT/ascat" "$PROJECT_OUT/sv/annotation" "$PROJECT_OUT/cohort" | awk '{print $4}')

if [ -f "${GROUP_TSV:-examples/group.example.tsv}" ]; then
  jid11=$(sbatch --dependency=afterok:"$jid10" slurm/group_compare.sbatch "${GROUP_TSV:-examples/group.example.tsv}" "$PROJECT_OUT/cohort" "$PROJECT_OUT/group_compare" | awk '{print $4}')
else
  jid11="$jid10"
fi

if [ -f "${HPV_BREAKPOINTS:-examples/hpv_breakpoints.example.tsv}" ]; then
  jid12=$(sbatch --dependency=afterok:"$jid11" slurm/hpv_link.sbatch "${HPV_BREAKPOINTS:-examples/hpv_breakpoints.example.tsv}" "$PROJECT_OUT/ascat" "$PROJECT_OUT/sv/annotation" "$PROJECT_OUT/hpv_link" | awk '{print $4}')
else
  jid12="$jid11"
fi

sbatch --dependency=afterok:"$jid12" slurm/final_report.sbatch "$PROJECT_OUT" "$PROJECT_OUT/final_report"
