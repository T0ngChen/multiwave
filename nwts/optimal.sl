#!/bin/bash
#SBATCH -J smalloptimal
#SBATCH -A uoa02593 # Project Account
#SBATCH --time=7:30:00 # Walltime
#SBATCH --mem=15GB # Memory per CPU
#SBATCH --cpus-per-task=5 # no. of core per task
#SBATCH --array=1-1000

module load JAGS/4.2.0-gimkl-2017a
module load R
# SLURM_JOBID
srun Rscript nwtsparany.r
