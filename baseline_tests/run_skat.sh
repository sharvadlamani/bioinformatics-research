#!/bin/bash

#SBATCH -N 1 # Ensure that all cores are on one machine
#SBATCH -p pe2,dev
#SBATCH --mem=10G
#SBATCH -t 0-23:00 # Runtime in D-HH:MM
#SBATCH -J SKAT # <-- name of job
#SBATCH --array=0-109# <-- number of jobs to run 
#SBATCH --output=/gpfs/commons/home/svadlamani/baseline_tests/outputs/stdout_%j.log               # Standard output and error log
#SBATCH --error=/gpfs/commons/home/svadlamani/baseline_tests/outputs/error_%j.log
#SBATCH --mail-type=FAIL                 # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=svadlamani@nygenome.org 

module load R/3.6.0

cells=("LHX2" "NEUN" "OLIG2" "peripheralPU1nuclei" "coding")
cell_index=$((SLURM_ARRAY_TASK_ID % ${#cells[@]}))
cell=${cells[$cell_index]}
echo $cell

chromosomes=(1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22)
chr_index=$((SLURM_ARRAY_TASK_ID % ${#chromosomes[@]}))
chr=${chromosomes[$chr_index]}
chr=$(((SLURM_ARRAY_TASK_ID / 5)%22 + 1)) 
echo $chr

iteration=$((SLURM_ARRAY_TASK_ID / 110 + 1))
echo $iteration

Rscript /gpfs/commons/home/svadlamani/baseline_tests/SKAT.R $chr $cell $iteration skat_all "ALL"
