#!/usr/bin/env bash
#SBATCH --partition=cpu1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=04:00:00
#SBATCH --job-name=sv_annotation
#SBATCH --output=logs/sv_annotation/%j.out
#SBATCH --error=logs/sv_annotation/%j.err

bash scripts/modules/sv_annotation.sh "$1" "$2" "$3" "$4"
