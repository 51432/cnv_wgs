#!/usr/bin/env bash
#SBATCH --partition=cpu1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=04:00:00
#SBATCH --job-name=group_compare
#SBATCH --output=logs/group_compare/%j.out
#SBATCH --error=logs/group_compare/%j.err

bash scripts/modules/group_compare.sh "$1" "$2" "$3"
