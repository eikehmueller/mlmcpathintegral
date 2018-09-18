#ifndef SAMPLER_HH
#define SAMPLER_HH SAMPLER_HH
#include <memory>
#include <iostream>
#include <iomanip>
#include <vector>
#include "path.hh"

/** @file sampler.hh
 * @brief Header file for sampler base class
 */

/** @class Sampler
 * @brief Base class for sampler
 *
 * Abstract base class for sampler from distribution over paths
 */
class Sampler {
public:
  /** @brief Create new instance
   *
   */
  Sampler() {
    reset_stats();
  }

  /** @brief Delete instance */
  virtual ~Sampler() {};

  /** @brief Draw a sample 
   *
   * returns a sample path \f$X\f$
   *
   * @param[out] x_path Path \f$X\f$ drawn from distribution
   */
  virtual void draw(std::vector<std::shared_ptr<Path>> x_path) = 0;

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
  
  
protected:
  /** @brief Collect statistics on acceptance probability and autocorrelation */
  mutable unsigned int n_accepted_samples;
  /** @brief Number of total samples */  
  mutable unsigned int n_total_samples;

};

#endif // SAMPLER_HH
