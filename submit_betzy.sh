#!/bin/bash
#SBATCH --job-name=demog1
## Project:
#SBATCH --account=nn9487k 
#SBATCH --nodes=4
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=128
#SBATCH --time=4-00:00:00
module restore system
#module load intel-compilers/2025.0.0
module load intel-compilers/2024.2.0
#module load mkl/2025.0
module load imkl/2024.2.0
module list
time ./MAIN > ProgramOutput.txt

