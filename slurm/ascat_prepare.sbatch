#!/usr/bin/env bash
#SBATCH --partition=cpu1
#SBATCH --cpus-per-task=8
#SBATCH --mem=48G
#SBATCH --time=12:00:00
#SBATCH --job-name=ascat_prepare
#SBATCH --output=logs/ascat_prepare/%A_%a.out
#SBATCH --error=logs/ascat_prepare/%A_%a.err

bash scripts/modules/ascat_prepare.sh "$1" "$2" "$3" "${SLURM_ARRAY_TASK_ID}"
