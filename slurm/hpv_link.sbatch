#!/usr/bin/env bash
#SBATCH --partition=cpu1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=04:00:00
#SBATCH --job-name=hpv_link
#SBATCH --output=logs/hpv_link/%A_%a.out
#SBATCH --error=logs/hpv_link/%A_%a.err

bash scripts/modules/hpv_link.sh "$1" "$2" "$3" "$4" "$5" "${SLURM_ARRAY_TASK_ID}"
