#!/bin/bash

#SBATCH --reservation=cpsc424
#SBATCH --cpus-per-task=20
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --time=4:00:00
#SBATCH --mem-per-cpu=5000
#SBATCH --job-name=ga-omp
#SBATCH --output=%x-%j.out

# Load Required Modules

module load intel

# set time format
export TIMEFORMAT=%2R

cp makefile_sort Makefile 
make clean
make 

echo ""
echo ""
echo "" OpenMP times
echo ""


for num_threads in 5
do
  export OMP_NUM_THREADS=$num_threads
  echo "Number of threads = " $OMP_NUM_THREADS
  time ./ga-omp-sort test
done
