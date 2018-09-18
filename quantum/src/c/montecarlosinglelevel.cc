#include "montecarlosinglelevel.hh"

/** @file montecarlosinglelevel.cc
 * @brief Implementation of montecarlosinglelevel.hh 
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
