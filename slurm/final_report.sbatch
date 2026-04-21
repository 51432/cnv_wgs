#!/usr/bin/env bash
#SBATCH --partition=cpu1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=02:00:00
#SBATCH --job-name=final_report
#SBATCH --output=logs/report/%j.out
#SBATCH --error=logs/report/%j.err

bash scripts/modules/final_report.sh "$1" "$2"
