#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <mpi.h>

#define FP double

int chromosomes;            // number of chromosomes
int genes;                  // number of genes per chromosome
FP beta;                    // variable for fitness function A
FP (*fitness)(int *);       // pointer to fitness function

// generate random integer in range [min, max]
int randint(int min, int max)
{
    return (rand() % (max - min + 1)) + min;
}

// generate random float in range [0, 1]
FP randf(void)
{
    return (FP) rand() / (FP) RAND_MAX;
}

// generate random integer in range [-genes, genes]
int random_gene(void)
{
    return randint(-genes, genes);
}

// get a pointer to the i-th individual in the population
int * get_individual(int * population, int i)
{
    return &population[i * genes];
}

// compute fitness according to fitness function A
FP fitness_A(int * individual)
{
    FP fit, sum, inner_sum;
    FP inner1, inner2;
    FP inner1_sq, inner2_sq;

    sum = 0.0;
    for (int k = 1; k <= genes; k++)
    {
        inner_sum = 0.0;
        for (int i = 1; i <= genes; i++)
        {
            inner1_sq = (FP) i;
            inner2_sq = (FP) individual[i - 1] / i;
            inner1 = 1;
            inner2 = 1;
            for (int j = 1; j <= k; j++)
            {
                inner1 *= inner1_sq;
                inner2 *= inner2_sq;
            }
            inner1 += (FP) beta;
            inner2 -= (FP) 1;
            inner_sum += inner1 * inner2;
        }
        sum += inner_sum * inner_sum;
    }
    fit = 1 - sum;

    return fit;
}

// compute fitness according to fitness function B
FP fitness_B(int * individual)
{
    FP fit, sum, inner_sum;
    FP inner1, inner2;
    FP inner1_sq, inner2_sq;

    sum = 0.0;
    for (int k = 1; k <= genes; k++)
    {
        inner_sum = 0.0;
        for (int i = 1; i <= genes; i++)
        {
            inner1_sq = (FP) individual[i - 1];
            inner2_sq = (FP) i;
            inner1 = 1;
            inner2 = 1;
            for (int j = 1; j <= k; j++)
            {
                inner1 *= inner1_sq;
                inner2 *= inner2_sq;
            }
            inner_sum += (inner1 - inner2);
        }
        sum += inner_sum * inner_sum;
    }
    fit = 1 - sum;

    return fit;
}

// print [size] members of [population]
void print_pop(int * population, int size)
{
    printf("[\n");
    for (int i = 0; i < size; i++) {
        printf("  [ ");
        for (int j = 0; j < genes; j++) {
            printf("%d, ", get_individual(population, i)[j]);
        }
        printf("\b\b ],\n");
    }
    printf("]\n");
}

/*
// calculate the average fitness of [size] members of [population]
FP average_fitness(int * population, int size)
{
    FP sum = 0.0;
    for (int i = 0; i < size; i++)
    {
        sum += fitness(get_individual(population, i));
    }
    return sum / (FP) size;
}
*/

// comparison function for qsort() fitness comparisons
int cmp_fitness(const void * chromosome1, const void * chromosome2)
{
    FP fit1 = fitness((int *) chromosome1);
    FP fit2 = fitness((int *) chromosome2);
    return (fit1 == fit2) ? 0 : (fit1 > fit2) ? -1 : 1;
}

// run one evolution iteration on the population
void evolve(int * population, int size, FP p_parents, FP p_crossover, FP p_mutate)
{
    // sort population by fitness
    qsort(population, size, genes * sizeof(int), cmp_fitness);

    // set number of parents
    int num_parents = p_parents * size;
    if ((size - num_parents) % 2 != 0) num_parents++;

    // create temporary variables
    int * parent1, * parent2, * child1, * child2;
    int parent1_gene, parent2_gene;
    int min, max, tmp_xi;

    // fill in the rest of the population with offspring of parents
    for (int chrom = num_parents; chrom < size; chrom += 2)
    {
        parent1 = get_individual(population, rand() % num_parents);
        parent2 = get_individual(population, rand() % num_parents);

        child1 = get_individual(population, chrom);
        child2 = get_individual(population, chrom + 1);

        for (int gene = 0; gene < genes; gene++)
        {
            parent1_gene = parent1[gene];
            parent2_gene = parent2[gene];

            // crossover the two parents in the children
            if (randf() < p_crossover)
            {
                child1[gene] = parent1_gene;
                child2[gene] = parent2_gene;
            }
            else
            {
                child1[gene] = parent2_gene;
                child2[gene] = parent1_gene;
            }

            // child1 random mutation
            if (randf() < p_mutate) child1[gene] = random_gene();

            // child2 parent-mixing mutation
            if (randf() < p_mutate)
            {
                if (parent1_gene != parent2_gene)
                {
                    min = (parent2_gene > parent1_gene) ? parent1_gene : parent2_gene;
                    max = (parent2_gene > parent1_gene) ? parent2_gene : parent1_gene;
                    child2[gene] = randint(min, max);
                }
                else child2[gene] = random_gene();
            }

            // child1 shifting mutation
            if (randf() < p_mutate)
            {
                tmp_xi = child1[gene];
                child1[gene] = child1[(gene - 1) % genes];
                child1[(gene - 1) % genes] = tmp_xi;
            }

            // child2 shifting mutation
            if (randf() < p_mutate)
            {
                tmp_xi = child2[gene];
                child2[gene] = child2[(gene - 1) % genes];
                child2[(gene - 1) % genes] = tmp_xi;
            }
        }
    }
}

/**
 * run_genetic_algorithm()
 * -----------------------
 * - runs a single instance of the genetic algorithm until termination
 * - assumes that the variables [genes], [chromosomes], and [fitness]
 *   are initialized before running
*/
void run_genetic_algorithm()
{
    // MPI variables
    int rank, procs;
    double wct1, wct2;

    // rank-0 variables
    int min, mod;
    int * population;
    int * proc_genes;
    int * proc_displs;
    int * best_of_the_best;

    // other variables
    int gen, max_gen;
    int island_genes;
    int island_chroms;
    int * island;
    int * best_island;
    int * best_overall;
    FP best_fitness;
    FP p_parents, p_crossover, p_mutate;

    // set random seed
    srand(12345);

    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &procs);

    if (rank == 0)
    {
        // allocate memory to hold the population
        population = (int *) calloc(chromosomes, genes * sizeof(int));

        // randomly initialize the population
        for (int chrom = 0; chrom < chromosomes; chrom++)
        {
            for (int gene = 0; gene < genes; gene++)
            {
                get_individual(population, chrom)[gene] = random_gene();
            }
        }

        // allocate memory to hold the displacement and size arrays
        proc_displs = (int *) calloc(procs, sizeof(int));
        proc_genes = (int *) calloc(procs, sizeof(int));

        // calculate displacement and size of each process
        min = chromosomes / procs;
        mod = chromosomes % procs;
        proc_displs[0] = 0;
        for (int proc = 0; proc < procs; proc++)
        {
            proc_genes[proc] = (min + ((proc >= mod) ? 0 : 1)) * genes;
            if (proc + 1 < procs) proc_displs[proc + 1] = proc_displs[proc] + proc_genes[proc];
        }

        // allocate array to hold fittest chromosome from each process
        best_overall = (int *) calloc(procs, genes * sizeof(int));
    }
    else
    {
        // allocate space to hold fittest chromosome from manager MPI_Bcast()
        best_overall = (int *) calloc(procs, genes * sizeof(int));
    }

    // send number of chromosomes to each process
    MPI_Scatter(proc_genes, 1, MPI_INT, &island_genes, 1, MPI_INT, 0, MPI_COMM_WORLD);

    // calculate number of chromosomes in each island
    island_chroms = island_genes / genes;

    // allocate enough space to hold required chromosomes
    island = (int *) calloc(island_chroms, genes * sizeof(int));

    // send the chromosomes themselves to each process
    MPI_Scatterv(population, proc_genes, proc_displs, MPI_INT,
                 island, island_genes, MPI_INT, 0, MPI_COMM_WORLD);

    // initialize
    gen = 0;
    max_gen = 10 * genes * genes * (2 * genes + 1);
    best_island = get_individual(island, 0);
    best_fitness = fitness(best_island);

    // evolution parameters
    p_parents = 0.3;
    p_crossover = 0.5;

    // adjust mutation probability based on number of chromosomes
    switch (chromosomes)
    {
        case 10:
            p_mutate = 0.5;
            break;
        case 50:
            p_mutate = 0.1;
            break;
        case 100:
            p_mutate = 0.2;
            break;
        case 150:
            p_mutate = 0.04;
            break;
        case 1000:
            p_mutate = 0.03;
        default:
            p_mutate = 0.05;
            break;
    }

    // synchronize before timing call
    MPI_Barrier(MPI_COMM_WORLD);
    if (rank == 0) wct1 = MPI_Wtime();

    // run evolution
    while (gen < max_gen)
    {
        // update generation
        gen++;

        // run a single evolution
        evolve(island, island_chroms, p_parents, p_crossover, p_mutate);

        // find fittest chromosome on the island
        best_island = get_individual(island, 0);

        // send fittest chromosome to manager process
        MPI_Gather(best_island, genes, MPI_INT, best_overall, genes, MPI_INT, 0, MPI_COMM_WORLD);

        if (rank == 0)
        {
            // sort fittest chromosomes from all regions
            qsort(best_overall, procs, genes * sizeof(int), cmp_fitness);
        }
        
        // broadcast fittest chromosome to all regions
        MPI_Bcast(best_overall, genes, MPI_INT, 0, MPI_COMM_WORLD);

        // copy overall fittest chromosome into island fittest chromosome
        best_overall = get_individual(best_overall, 0);
        for (int gene = 0; gene < genes; gene++) best_island[gene] = best_overall[gene];

        // update best fitness
        best_fitness = fitness(best_island);

        // check exit condition
        if (best_fitness == 1.0) break;
    }

    // collect the population back to the root process
    MPI_Gatherv(island, island_genes, MPI_INT,
                population, proc_genes, proc_displs, MPI_INT, 0, MPI_COMM_WORLD);

    // synchronize before manager-only section
    MPI_Barrier(MPI_COMM_WORLD);

    if (rank == 0)
    {
        qsort(population, chromosomes, genes * sizeof(int), cmp_fitness);

        // get end time
        wct2 = MPI_Wtime();

        // find best individual
        best_of_the_best = get_individual(population, 0);
        best_fitness = fitness(best_of_the_best);

        // print overall statistics
        printf("==========================================\n");
        printf("Number of MPI Processes:      %12d\n", procs);
        printf("Computation Time:             %12.8f\n", wct2 - wct1);
        printf("==========================================\n");
        printf("Population Size:              %12d\n", chromosomes);
        printf("Number of Generations:        %12d\n", gen);
        printf("==========================================\n");
        printf("Best Individual Fitness:   %+12.8e\n", best_fitness);
        printf("Best Individual: [");
        for (int gene = 0; gene < genes; gene++)
        {
            printf("%d", best_of_the_best[gene]);
            if (gene + 1 != genes) printf(", ");
            else printf("]\n");
        }
        printf("==========================================\n");

        // free manager-only allocations
        free(population);
        free(proc_genes);
        free(proc_displs);
    }
    // free allocations in MPI processes
    free(best_overall);
    free(island);
}


int main(int argc, char ** argv)
{
    if (argc != 5)
    {
        fprintf(stderr, "Usage: ./ga-island <chromosomes> <genes> <beta> <fitness function>\n");
        exit(1);
    }

    chromosomes = atoi(argv[1]);
    genes = atoi(argv[2]);
    beta = atof(argv[3]);
    fitness = (strcmp(argv[4], "A") == 0) ? fitness_A : fitness_B;

    MPI_Init(&argc, &argv);

    run_genetic_algorithm();

    MPI_Finalize();

    return 0;
}
