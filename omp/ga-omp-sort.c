#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <string.h> // memmove
#include <math.h>
#include <omp.h>
#include "timing.h"

#define FP double

int population_size;
int n;
FP beta;
FP (*fitness)(int *);

// returns a random point in range [-n, n]
int random_xi(void) {
    return (rand() % (2*n+1)) - n;   
}

int randint(int min, int max) {
    return (rand() % (max - min + 1)) + min;
}

// return a pointer to the i-th individual
int * get_individual(int * pop, int i) {
    return &pop[i*n];
}

void swap_individuals(int * pop, int ind1_i, int ind2_i) {
  int temp;
  int * ind1 = get_individual(pop, ind1_i);
  int * ind2 = get_individual(pop, ind2_i);
  for (int j=0; j < n; j++) {
      temp = ind1[j];
      ind1[j] = ind2[j];
      ind2[j] = temp;
  }
}

int partition(int * pop, int low, int high) {
  FP pivot = fitness(get_individual(pop, high));
  int i = low-1; 
  for (int j = low; j <= high-1; j++) {
      if (fitness(get_individual(pop, j)) > pivot) {
          i++; 
          swap_individuals(pop, i, j);
      }
  }
  swap_individuals(pop, i+1, high);
  return i+1;
}

void quicksort(int * pop, int low, int high) {
  if (low < high) {
    int pivot_i = partition(pop, low, high);

    #pragma omp task
    quicksort(pop, low, pivot_i-1);

    #pragma omp task 
    quicksort(pop, pivot_i+1, high);
  }
}

// calculate the individual fitness of the individual
FP fitness_A(int * individual)
{
    FP fit, sum, inner_sum;
    FP inner1, inner2;
    FP inner1_sq, inner2_sq;

    sum = 0.0;
    for (int k = 1; k <= n; k++)
    {
        inner_sum = 0.0;
        for (int i = 1; i <= n; i++)
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

// calculate the individual fitness of the individual
FP fitness_B(int * individual)
{
    FP fit, sum, inner_sum;
    FP inner1, inner2;
    FP inner1_sq, inner2_sq;

    sum = 0.0;
    for (int k = 1; k <= n; k++)
    {
        inner_sum = 0.0;
        for (int i = 1; i <= n; i++)
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

void print_pop(int * pop) {
    printf("[\n");
    for (int i = 0; i < population_size; i++) {
        printf("  [ ");
        for (int j = 0; j < n; j++) {
            printf("%d, ", get_individual(pop, i)[j]);
        }
        printf("\b\b ],\n");
    }
    printf("]\n");
}

// calculate the average fitness of the population
FP average_fitness(int * pop) {
    FP sum = 0.0;
    int i;
    for (i = 0; i < population_size; i++) {
        sum += fitness(get_individual(pop, i));
    }
    return sum / (FP) population_size;
}

int cmp_fitness (const void * a, const void * b) {
   return (int) (fitness((int *) b) - fitness((int *) a));
}
        
// returns a random float from 0 to 1
FP randf(void) {
    return (FP) rand() / (FP) RAND_MAX;
}

void evolve_pop(int * pop, FP p_mutate, FP p_elite, FP p_crossover, FP p_parents) {
    int j, parent1, parent2, parent1_xi, parent2_xi, old_xi;
    int * child1, * child2;

    //int num_elite = (p_elite * population_size >= 1) ? p_elite*population_size : (p_elite > 0);

    int num_parents = p_parents * population_size;
    if ((population_size - num_parents) % 2 != 0) num_parents++;

    int next_row = num_parents;

    // sort pop by fitness
    #pragma omp parallel
    {
      #pragma omp single
      quicksort(pop, 0, population_size-1);
    }

    // the rest of the popluation will consist of a breed between 
    // two random parents
    for (; next_row < population_size; next_row += 2) {
        parent1 = rand() % num_parents;
        while ((parent2 = rand() % num_parents) == parent1)
            ; // wait till we get a different parent
            
        child1 = get_individual(pop, next_row);
        child2 = get_individual(pop, next_row+1);

        for (j = 0; j < n; j++) {
            // create the children
            parent1_xi = get_individual(pop, parent1)[j];
            parent2_xi = get_individual(pop, parent2)[j];
            if (randf() < p_crossover) {
                child1[j] = parent1_xi;
                child2[j] = parent2_xi;
            } else {
                child1[j] = parent2_xi;
                child2[j] = parent1_xi;
            }

            // mutate the child
            if (randf() < p_mutate) { // mutate child1 randomly
                child1[j] = random_xi();
            }
            if (randf() < p_mutate) { // mutate child2 in between parents
                if (parent1_xi != parent2_xi)
                    child2[j] = randint(((parent2_xi > parent1_xi) ? parent1_xi : parent2_xi),
                            ((parent2_xi > parent1_xi) ? parent2_xi : parent1_xi));
                else
                    child2[j] = random_xi();
            }

            if (randf() < p_mutate) { // randomly shift digits
                old_xi = child1[j];
                child1[j] = child1[(j - 1) % n];
                child1[(j - 1) % n] = old_xi;
            }
            if (randf() < p_mutate) { // mutate child2 randomly
                old_xi = child2[j];
                child2[j] = child2[(j - 1) % n];
                child2[(j - 1) % n] = old_xi;
            }
        }
    }
}

FP run_ga(int num_trials) {
    int i;
    int max_iters = 10*n*n*(2*n+1);
    int * population = NULL;
    int * best_individual = NULL;
    FP indiv_best_fitness = 0.;
    FP start_time = 0.;
    FP end_time = 0.;
    FP total_time = 0.;
    long total_gens = 0;

    srand(12345);

    for (int trial = 0; trial < num_trials; trial++) {
        // allocate array of individuals
        population = (int *) malloc(n * population_size * sizeof(int));

        // initialize individuals randomly
        for (i = 0; i < population_size; i++)
        {
            for (int xi = 0; xi < n; xi++)
            {
                if (fitness == fitness_B && beta == 1 && i < population_size/20)
                    get_individual(population, i)[xi] = xi;
                else 
                    get_individual(population, i)[xi] = random_xi();
            }
        }

        FP best_fitness = fitness(get_individual(population, 0));      // current optimal fitness

        FP p_mutate = 0.05;
        FP p_elite = 0.01;
        FP p_crossover = 0.5;
        FP p_parents = 0.3;
        // adjust mutation probability based on population size
        switch (population_size) {
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
                break;
        }
                                
        timing(&start_time, NULL);
        int gen = 0;            // current generation
        while (gen < max_iters)
        {
            // update current generation
            gen += 1;
            
            // evolve population
            evolve_pop(population, p_mutate, p_elite, p_crossover, p_parents);
            //print_pop(population); 

            best_fitness = fitness(get_individual(population, 0));
            if (best_fitness == 1.0) break; // reached them maximum
        }

        // sort pop by fitness one last time
        #pragma omp parallel
        {
          #pragma omp single
          quicksort(population, 0, population_size-1);
        }

        best_individual = get_individual(population, 0);
        indiv_best_fitness = fitness(best_individual);

        if (trial != num_trials - 1) {
            free(population);
        }
        timing(&end_time, NULL);
        total_time += (end_time - start_time);
        total_gens += gen;
    }
        
    // print overall statistics
    printf("========================================\n");
    printf("Population Size:            %12d\n", population_size);
    printf("Number of Trials:           %12d\n", num_trials);
    printf("Avg Number of Generations:  %12f\n", (FP) total_gens / (FP) num_trials);
    printf("Avg Time (ms):              %12f\n", total_time / (FP) num_trials);
    printf("Best Individual Fitness: %+12.8e\n", indiv_best_fitness);
    printf("Best Individual:         [");
    for (i = 0; i < n-1; i++) {
        printf("%d, ", best_individual[i]);
    }
    printf("%d]\n", best_individual[n-1]);
    printf("========================================\n");

    free(population);

    return indiv_best_fitness;    
}

#define NUM_TESTS 4
#define NUM_POP_SIZES 5

int main(int argc, char ** argv)
{
    char fitness_letter;
    if (argc == 5) {
        population_size = atoi(argv[1]);
        n = atoi(argv[2]);
        beta = atof(argv[3]);
        fitness_letter = argv[4][0];
        
        fitness = (fitness_letter == 'A') ? fitness_A : fitness_B;
        
        run_ga(1);
    } else if (argc == 2) {
        int A_ns[NUM_TESTS] = {4, 4, 10, 10};
        FP A_betas[NUM_TESTS] = {50, 0.5, 1.0e9, 1.0e7};
        int B_ns[NUM_TESTS] = {5, 5, 10, 10};
        FP B_betas[NUM_TESTS] = {0, 1, 0, 1};
        int pop_sizes[NUM_POP_SIZES] = { 10, 50, 100, 150, 1000 };
        int pop_i;
        for (pop_i = 0; pop_i < NUM_POP_SIZES; pop_i++) {
            population_size = pop_sizes[pop_i];
            for (int trial = 0; trial < NUM_TESTS-2; trial++) {
                fitness = fitness_A;
                n = A_ns[trial]; beta = A_betas[trial];
                printf("\n\nA with pop size: %d, n: %d, beta: %.1f\n", population_size, n, beta);
                run_ga(3);

                fitness = fitness_B;
                n = B_ns[trial]; beta = B_betas[trial];
                printf("\n\nB with pop size: %d, n: %d, rand: %.1f\n", population_size, n, beta);
                run_ga(3);
            }
        }
    } else {
        fprintf(stderr, "Usage: ./ga <population size> <n> <beta> <fitness func letter>\n");
        return 1;
    } 

    return 0;
}
