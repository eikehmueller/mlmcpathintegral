#include "montecarlo.hh"
/** @brief Implementation of montecarlo.hh */

/** Calculate Monte Carlo estimate with single level method */
std::pair<double,double> MonteCarloSingleLevel::evaluate() {
  std::vector<Path*> x_path(1);
  x_path[0] = new Path(sampler.getM_lat());

  // Burn-in phase
  for (unsigned int k=0;k<n_burnin;++k) {
    sampler.draw(x_path);
  }

  // Sampling phase
  double S = 0.0;
  double S2 = 0.0;
  for (unsigned int k=0;k<n_samples;++k) {
    sampler.draw(x_path);
    double tmp=qoi.evaluate(x_path[0]);
    S += tmp;
    S2 += tmp*tmp;
  }
  double mean = S/n_samples;
  double error = sqrt(1./(n_samples*(n_samples-1.))*(S2-S*S/n_samples));
  return std::make_pair(mean,error);
}

/** Calculate mean and variance of difference in QoI */
std::pair<double,double> MonteCarloTwoLevel::evaluate_difference() {
  std::vector<Path*> x_path(2);
  x_path[0] = new Path(fine_action.getM_lat()); // fine path
  x_path[1] = new Path(coarse_action.getM_lat()); // coarse path

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
