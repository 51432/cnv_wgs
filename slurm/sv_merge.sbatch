#!/usr/bin/env bash
#SBATCH --partition=cpu1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=02:00:00
#SBATCH --job-name=sv_merge
#SBATCH --output=logs/sv_merge/%A_%a.out
#SBATCH --error=logs/sv_merge/%A_%a.err

bash scripts/modules/sv_merge.sh "$1" "$2" "$3" "${SLURM_ARRAY_TASK_ID}"
