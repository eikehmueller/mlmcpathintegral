#include "mpi_random.hh"
/** @file Implementation of mpi_random.hh */

/* Seed random number generator */
void parallel_mt19937_64::seed(result_type value) {
  int n_rank = mpi_comm_size(); // Number of processors
  unsigned int* master_seed_list;
  unsigned int local_seed; // Local seed
  typedef std::linear_congruential_engine<unsigned int, 48271, 0, 2147483647> SeedEngine;
  // Generate seeds on master processor
  if (mpi_master()) {
    SeedEngine seed_engine;
    seed_engine.seed(value);
    std::set<SeedEngine::result_type> seed_set;
    seed_set.insert(value);
    while (seed_set.size() < n_rank) {
      seed_set.insert(seed_engine());
    }
    master_seed_list = (unsigned int*) malloc(n_rank*sizeof(unsigned int));
    std::copy(seed_set.begin(),seed_set.end(),master_seed_list);
  }
  // Distribute seeds to all processes
  mpi_scatter(master_seed_list,local_seed);
  if (mpi_master()) {
    free(master_seed_list);
  }
  base_engine.seed(local_seed);
}
