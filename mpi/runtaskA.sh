#!/bin/bash
#SBATCH --reservation=cpsc424
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --cpus-per-task=1
#SBATCH --tasks-per-node=4
#SBATCH --mem-per-cpu=5G
#SBATCH --time=20:00
#SBATCH --job-name=GA-MPI-islandA
#SBATCH --output=%x-%j.out

# Load necessary module files
module load iomkl
module list

# Print initial working directory
echo
echo -n "Working Directory: "
pwd

echo
echo "Making ga-island"
make clean
make ga-island

# Print the node list
echo
echo
echo "Node List:"
echo $SLURM_NODELIST
echo "ntasks-per-node = " $SLURM_NTASKS_PER_NODE

# Run the program on all required inputs
for chromosomes in 50 100 150 1000
do
    for genes in 4 8
    do
        for beta in 50 0.5
        do
            echo "chromosomes =" $chromosomes
            echo "genes =" $genes
            echo "beta =" $beta
            time mpirun -n 4 ga-island $chromosomes $genes $beta A
	    echo
	    echo
        done
    done
done
