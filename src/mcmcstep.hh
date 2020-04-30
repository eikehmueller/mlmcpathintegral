#ifndef MCMCSTEP_HH
#define MCMCSTEP_HH MCMCSTEP_HH
#include <iostream>
#include <iomanip>
#include <vector>
#include "mpi_wrapper.hh"

/** @file mcmcstep.hh
 * @brief Header file for MCMCStep class
 */

/** @class MCMCStep
 * @brief Base class for MCMC type samplers
 *
 * Base class for samples which transition from one state to the
 * next in an accept/reject step. This class mainly provides functionality
 * for recording the acceptance probability
 */
class MCMCStep {
public:
  /** @brief Create new instance
   *
   */
  MCMCStep() : accept(false) {
    reset_stats();
  }

  /** @brief reset statistics
   * 
   * Reset all sampling statistics
   */
  void reset_stats() {
    n_total_samples = 0;
    n_accepted_samples = 0;
  }

  /** @brief Return acceptance probability */
  double p_accept() { return n_accepted_samples/(1.*n_total_samples); }

  /** @brief Show statistics 
   *
   * Print out statistics 
   */
  void show_stats();

  /** @brief has last sample been accepted? */
  bool accepted() const { return accept; }
  
protected:
  /** @brief Collect statistics on acceptance probability and autocorrelation */
  mutable unsigned int n_accepted_samples;
  /** @brief Number of total samples */  
  mutable unsigned int n_total_samples;
  /** @brief Has last sample been accepted? */
  mutable bool accept;
};

#endif // MCMCSTEP_HH
