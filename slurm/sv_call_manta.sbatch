#!/usr/bin/env bash
#SBATCH --partition=cpu1
#SBATCH --cpus-per-task=8
#SBATCH --mem=48G
#SBATCH --time=10:00:00
#SBATCH --job-name=sv_manta
#SBATCH --output=logs/manta/%A_%a.out
#SBATCH --error=logs/manta/%A_%a.err

bash scripts/modules/sv_call_manta.sh "$1" "$2" "$3" "${SLURM_ARRAY_TASK_ID}"
