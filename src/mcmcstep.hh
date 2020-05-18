#ifndef MCMCSTEP_HH
#define MCMCSTEP_HH MCMCSTEP_HH
#include <iostream>
#include <iomanip>
#include <vector>
#include "mpi_wrapper.hh"
#include "path.hh"

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
  MCMCStep() :
    accept(false),
    copy_if_rejected(false) {
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
  virtual void show_stats();

  /** @brief has last sample been accepted? */
  bool accepted() const { return accept; }
  
  /** @brief Set current state to particular value
   *
   * @param[in] x_path
   */
  virtual void set_state(std::shared_ptr<Path> x_path)=0;
  
  /** Return cost per sample */
  virtual double cost_per_sample() {
    mpi_parallel::cerr << " ERROR: Cost per sample not defined for class " << typeid(*this).name() << std::endl;
    mpi_exit(EXIT_FAILURE);
    return 0.0;
  }
  
protected:
  /** @brief Collect statistics on acceptance probability and autocorrelation */
  mutable unsigned int n_accepted_samples;
  /** @brief Number of total samples */  
  mutable unsigned int n_total_samples;
  /** @brief Has last sample been accepted? */
  mutable bool accept;
  /** @brief Copy path even if it has been rejected */
  mutable bool copy_if_rejected;
};

#endif // MCMCSTEP_HH
