#!/usr/bin/env bash
#SBATCH --partition=cpu1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=02:00:00
#SBATCH --job-name=input_check
#SBATCH --output=logs/input_check/%j.out
#SBATCH --error=logs/input_check/%j.err

bash scripts/modules/input_check.sh "$1" "$2" "$3"
