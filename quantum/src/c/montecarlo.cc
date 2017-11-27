#include "montecarlo.hh"
/** @brief Implementation of montecarlo.hh */
std::pair<double,double> MonteCarloSingleLevel::evaluate() {
  std::vector<double*> x_path(1);
  x_path[0] = new double[sampler.getM_lat()];

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
  delete[] x_path[0];
  return std::make_pair(mean,error);
}
