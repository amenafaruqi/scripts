#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=48
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=3700
#SBATCH --time=48:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=amena.faruqi@warwick.ac.uk
#SBATCH --job-name="ODMS"

module purge
module restore
module load intel


MY_NUM_THREADS=$SLURM_NTASKS_PER_NODE

export OMP_NUM_THREADS=$MY_NUM_THREADS

srun -n 1 -c $MY_NUM_THREADS python3 generate_fits.py
