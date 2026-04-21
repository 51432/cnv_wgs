#!/usr/bin/env bash
#SBATCH --partition=cpu1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=12:00:00
#SBATCH --job-name=ascat_run
#SBATCH --output=logs/ascat/%A_%a.out
#SBATCH --error=logs/ascat/%A_%a.err

bash scripts/modules/ascat_run.sh "$1" "$2" "$3" "$4" "${SLURM_ARRAY_TASK_ID}"
