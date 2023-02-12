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

cp makefile_fit Makefile
make clean
make

echo ""
echo ""
echo "" OpenMP times
echo ""


for schedule in "static,1" "dynamic" "guided"
do
    # OpenMP Parameter Settings
    export OMP_SCHEDULE=$schedule
    echo "OMP_SCHEDULE = " $OMP_SCHEDULE

     for num_threads in 4 5
     do
         export OMP_NUM_THREADS=$num_threads
         echo "Number of threads = " $OMP_NUM_THREADS
         time ./ga-omp-fit test
     done
 done
