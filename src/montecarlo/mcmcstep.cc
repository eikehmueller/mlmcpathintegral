#include "mcmcstep.hh"

/** @file mcmcstep.cc
 * @brief Implementation of mcmcstep.hh
 */

/** Show statistics */
void MCMCStep::show_stats() {
    mpi_parallel::cout << std::setprecision(3) << std::fixed;
    mpi_parallel::cout << "  acceptance probability  p = " << p_accept() << std::endl;
    mpi_parallel::cout << "  rejection probability 1-p = " << 1.-p_accept() << std::endl;
}
