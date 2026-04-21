#!/usr/bin/env bash
#SBATCH --partition=cpu1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=04:00:00
#SBATCH --job-name=cohort_summary
#SBATCH --output=logs/cohort_summary/%j.out
#SBATCH --error=logs/cohort_summary/%j.err

bash scripts/modules/cohort_summary.sh "$1" "$2" "$3"
