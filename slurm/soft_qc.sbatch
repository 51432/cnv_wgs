#!/usr/bin/env bash
#SBATCH --partition=cpu1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=04:00:00
#SBATCH --job-name=soft_qc
#SBATCH --output=logs/soft_qc/%A_%a.out
#SBATCH --error=logs/soft_qc/%A_%a.err

bash scripts/modules/soft_qc.sh "$1" "$2" "${SLURM_ARRAY_TASK_ID}"
