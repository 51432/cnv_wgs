#!/usr/bin/env bash
#SBATCH --partition=cpu1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=02:00:00
#SBATCH --job-name=sv_merge
#SBATCH --output=logs/sv_merge/%j.out
#SBATCH --error=logs/sv_merge/%j.err

bash scripts/modules/sv_merge.sh "$1" "$2" "$3"
