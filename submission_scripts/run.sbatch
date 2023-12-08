#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=3700
#SBATCH --time=00:10:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=amena.faruqi@warwick.ac.uk
#SBATCH --job-name=test

#ulimit -s unlimited
MY_NUM_THREADS=$SLURM_CPUS_PER_TASK

export OMP_NUM_THREADS=$MY_NUM_THREADS
#export OMP_DYNAMIC=FALSE
#export OMP_MAX_ACTIVE_LEVELS=1
LAST_OUTPUT=$(ls -Art gasdens* | grep -Eo '[0-9]+' | tail -n 1)
echo ${LAST_OUTPUT}

srun $FARGODIR/fargo3d -S $LAST_OUTPUT param.par > output

