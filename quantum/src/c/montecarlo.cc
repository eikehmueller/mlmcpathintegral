#include "montecarlo.hh"
#include "config.h"
#include <sstream>
#include <iomanip>

/** @file montecarlo.cc
 * @brief Implementation of montecarlo.hh 
 */

/** Calculate Monte Carlo estimate with single level method */
std::pair<double,double> MonteCarloSingleLevel::evaluate() {
  std::vector<Path*> x_path(1);
  x_path[0] = new Path(action.getM_lat(),
                       action.getT_final());

  // Burn-in phase
  for (unsigned int k=0;k<n_burnin;++k) {
    sampler.draw(x_path);
  }

  // Sampling phase
  double S = 0.0;
  double S2 = 0.0;
  for (unsigned int k=0;k<n_samples;++k) {
    sampler.draw(x_path);
    /* Save (some) paths to disk? Edit file config.h */
#ifdef SAVE_PATHS
    if ( (SAVE_FIRST_PATH<=k) and (k<=SAVE_LAST_PATH) ) {
      std::stringstream filename;      
      filename << "path_";
      filename << std::setw(8) << std::setfill('0');
      filename << k << ".dat";
      x_path[0]->save_to_disk(filename.str());
    }
#endif
    double tmp=qoi.evaluate(x_path[0]);
    S += tmp;
    S2 += tmp*tmp;
  }
  double mean = S/n_samples;
  double variance = (S2-S*S/n_samples)/(n_samples-1.);
  return std::make_pair(mean,variance);
}

/** Calculate mean and variance of difference in QoI */
std::pair<double,double> MonteCarloTwoLevel::evaluate_difference() {
  std::vector<Path*> x_path(2);
  x_path[0] = new Path(fine_action.getM_lat(),
                       fine_action.getT_final()); // fine path
  x_path[1] = new Path(coarse_action.getM_lat(),
                       fine_action.getT_final()); // coarse path

  // Burn-in phase
  for (unsigned int k=0;k<n_burnin;++k) {
    twolevel_sampler.draw(x_path);
  }

  // Sampling phase
  double S = 0.0;
  double S2 = 0.0;
  for (unsigned int k=0;k<n_samples;++k) {
    twolevel_sampler.draw(x_path);
    double tmp=qoi.evaluate(x_path[0])-qoi.evaluate(x_path[1]);
    S += tmp;
    S2 += tmp*tmp;
  }
  double mean = S/n_samples;
  double variance = 1./(n_samples-1.)*(S2-S*S/n_samples);
  return std::make_pair(mean,variance);
}
