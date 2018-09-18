#include "montecarlo.hh"
#include "config.h"
#include <sstream>
#include <iomanip>

/** @file montecarlo.cc
 * @brief Implementation of montecarlo.hh 
 */

/** Calculate Monte Carlo estimate with single level method */
void MonteCarloSingleLevel::evaluate(Statistics& stats) {
  std::shared_ptr<Path> x_path =
    std::make_shared<Path>(action->getM_lat(),
                           action->getT_final());

  // Burn-in phase
  for (unsigned int k=0;k<n_burnin;++k) {
    sampler->draw(x_path);
  }
  stats.reset();
  for (unsigned int k=0;k<n_samples;++k) {
    sampler->draw(x_path);
    /* Save (some) paths to disk? Edit file config.h */
#ifdef SAVE_PATHS
    if ( (SAVE_FIRST_PATH<=k) and (k<=SAVE_LAST_PATH) ) {
      std::stringstream filename;      
      filename << "path_";
      filename << std::setw(8) << std::setfill('0');
      filename << k << ".dat";
      x_path->save_to_disk(filename.str());
    }
#endif
    stats.record_sample(qoi->evaluate(x_path));
  }
}

/** Calculate mean and variance of difference in QoI */
void MonteCarloTwoLevel::evaluate_difference(Statistics& stats_fine,
                                             Statistics& stats_coarse,
                                             Statistics& stats_diff) {
  std::shared_ptr<Path> x_path =
    std::make_shared<Path>(fine_action->getM_lat(),
                           fine_action->getT_final());
  std::shared_ptr<Path> x_coarse_path =
    std::make_shared<Path>(coarse_action->getM_lat(),
                           coarse_action->getT_final());

  // Burn-in phase
  for (unsigned int k=0;k<n_burnin;++k) {
    coarse_sampler->draw(x_coarse_path);
    twolevel_step.draw(x_coarse_path,x_path);
  }

  // Sampling phase
  stats_diff.reset();
  for (unsigned int k=0;k<n_samples;++k) {
    coarse_sampler->draw(x_coarse_path);
    twolevel_step.draw(x_coarse_path,x_path);
    double qoi_fine = qoi->evaluate(x_path);
    double qoi_coarse = qoi->evaluate(x_coarse_path);
    stats_fine.record_sample(qoi_fine);
    stats_coarse.record_sample(qoi_coarse);
    stats_diff.record_sample(qoi_fine-qoi_coarse);
  }
}
