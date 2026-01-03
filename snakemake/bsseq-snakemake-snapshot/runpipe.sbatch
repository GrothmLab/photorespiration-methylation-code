#!/bin/bash

################ sbatch parameters ################

#SBATCH -o %x_%j.out
#SBATCH -e %x_%j.err
#SBATCH -J snake_WGBS
#SBATCH -q cpu_normal #replace 
#SBATCH -p cpu_p #replace
#SBATCH --mem=500M

### task:

exec snakemake --profile my_snakemake_slurm_profile --rerun-triggers mtime "$@"
