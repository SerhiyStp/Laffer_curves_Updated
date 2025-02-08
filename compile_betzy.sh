#!/bin/bash
module restore system
#module load intel-compilers/2025.0.0
module load intel-compilers/2024.2.0
#module load mkl/2025.0
module load imkl/2024.2.0
module list
make
