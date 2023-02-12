#!/bin/bash
#SBATCH --reservation=cpsc424
#SBATCH --nodes=1
#SBATCH --ntasks=5
#SBATCH --cpus-per-task=1
#SBATCH --tasks-per-node=10
#SBATCH --mem-per-cpu=5G
#SBATCH --time=20:00
#SBATCH --job-name=GA-MPI-islandB
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
    for genes in 5 10
    do
        echo "chromosomes =" $chromosomes
        echo "genes =" $genes
        time mpirun -n 5 ga-island $chromosomes $genes 0 B
	echo
	echo
    done
    if [ $chromosomes -eq 1000 ]
    then
	for genes in 15 20 25
	do
	    echo "chromosomes =" $chromosomes
            echo "genes =" $genes
            time mpirun -n 5 ga-island $chromosomes $genes 0 B
            echo
            echo
	done
    fi
done
