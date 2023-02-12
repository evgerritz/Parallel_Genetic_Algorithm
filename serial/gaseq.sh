#!/bin/bash

#SBATCH --reservation=cpsc424
#SBATCH --cpus-per-task=20
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --time=4:00:00
#SBATCH --mem-per-cpu=5000
#SBATCH --job-name=gaseq
#SBATCH --output=%x-%j.out

# This script loads module files, builds programs, and then runs the programs for Assignment 2.
#    (Note it is only set up for a couple of the tasks, but you can modify it or enlarge it 
#     to encompass the other tasks, as well.)

# The script uses the Makefile and drand.c files in /home/cpsc424_ahs3/shared/a2. 

# To run the script, submit it to Slurm using: sbatch slurmrun.sh. Note that it requests a full node (20 cores).
# During program development, you may want a simplified version of the script (such as commenting out 
# or removing unnecessary lines. (You may prefer to use an interactive remote desktop session during development.)

# The script will produce an output file named something like "Mandelbrot1-2-49530642.out"
#    where 49530642 is the Slurm job number.

# Load Required Modules

module load intel

# set time format
export TIMEFORMAT=%2R

make clean
make

# Task 1

echo ""
echo ""
echo "Serial version"

time ./ga test

