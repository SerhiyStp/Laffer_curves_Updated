#!/bin/bash
#SBATCH --job-name=demog1
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=192
#SBATCH --time=2-10:00:00
module restore system
module load intel-compilers/2025.0.0
module load mkl/2025.0
module list
time ./MAIN > ProgramOutput.txt

