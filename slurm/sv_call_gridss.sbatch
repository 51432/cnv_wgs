#!/usr/bin/env bash
#SBATCH --partition=cpu2
#SBATCH --cpus-per-task=24
#SBATCH --mem=192G
#SBATCH --time=24:00:00
#SBATCH --job-name=sv_gridss
#SBATCH --output=logs/gridss/%A_%a.out
#SBATCH --error=logs/gridss/%A_%a.err

bash scripts/modules/sv_call_gridss.sh "$1" "$2" "$3" "${SLURM_ARRAY_TASK_ID}"
