#!/usr/bin/env bash
#SBATCH --partition=cpu1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=04:00:00
#SBATCH --job-name=sv_postfilter
#SBATCH --output=logs/sv_postfilter/%A_%a.out
#SBATCH --error=logs/sv_postfilter/%A_%a.err

bash scripts/modules/sv_postfilter.sh "$1" "$2" "$3" "$4" "${SLURM_ARRAY_TASK_ID}"
