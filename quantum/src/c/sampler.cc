#include "sampler.hh"

/** @file sampler.cc
 * @brief Implementation of sampler.hh
 */

/** Show statistics */
void Sampler::show_stats() {
  if (record_stats) {
    std::cout << "  acceptance probability  p = " << p_accept() << std::endl;
    std::cout << "  rejection probability 1-p = " << 1.-p_accept() << std::endl;
  } else {
    std::cout << "  No sampling statistics available." << std::endl;
  }
}

