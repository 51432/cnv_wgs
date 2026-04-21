#!/usr/bin/env bash
#SBATCH --partition={{partition}}
#SBATCH --cpus-per-task={{cpus_per_task}}
#SBATCH --mem={{mem}}
#SBATCH --time={{time}}
#SBATCH --job-name={{job_name}}
#SBATCH --array=1-{{array_size}}%{{array_parallel}}
#SBATCH --output={{stdout}}
#SBATCH --error={{stderr}}

{{command}}
