#!/usr/bin/env bash
#SBATCH --partition={{partition}}
#SBATCH --cpus-per-task={{cpus_per_task}}
#SBATCH --mem={{mem}}
#SBATCH --time={{time}}
#SBATCH --job-name={{job_name}}
#SBATCH --output={{stdout}}
#SBATCH --error={{stderr}}

{{command}}
