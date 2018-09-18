#include "mcmcstep.hh"

/** @file mcmcstep.cc
 * @brief Implementation of mcmcstep.hh
 */

/** Show statistics */
void MCMCStep::show_stats() {
  std::cout << std::setprecision(3) << std::fixed;
  std::cout << "  acceptance probability  p = " << p_accept() << std::endl;
  std::cout << "  rejection probability 1-p = " << 1.-p_accept() << std::endl;
}

