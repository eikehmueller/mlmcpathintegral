#include "sampler.hh"

/** @file sampler.cc
 * @brief Implementation of sampler.hh
 */

/** Show statistics */
void Sampler::show_stats() {
  std::cout << std::setprecision(3) << std::fixed;
  std::cout << "  acceptance probability  p = " << p_accept() << std::endl;
  std::cout << "  rejection probability 1-p = " << 1.-p_accept() << std::endl;
}

