#!/usr/bin/env bash
#SBATCH --partition=cpu1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=04:00:00
#SBATCH --job-name=sv_annotation
#SBATCH --output=logs/sv_annotation/%A_%a.out
#SBATCH --error=logs/sv_annotation/%A_%a.err

bash scripts/modules/sv_annotation.sh "$1" "$2" "$3" "$4" "$5" "${SLURM_ARRAY_TASK_ID}"
