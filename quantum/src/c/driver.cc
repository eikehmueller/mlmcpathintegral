#include <iostream>
#include <utility>
#include <memory>
#include "harmonicoscillatoraction.hh"
#include "quantityofinterest.hh"
#include "montecarlo.hh"
#include "parameters.hh"
#include "renormalisation.hh"
#include "twolevelmetropolissampler.hh"
#include "hmcsampler.hh"

/** @file driver.cc
 * @brief File with main program
 */

/** Main program */
int main(int argc, char* argv[]) {
  std::cout << "*** Path integral multilevel MCMC ***" << std::endl;  
  Parameters param("parameters.in");
  param.show();
  QoIXsquared qoi(param.M_lat);
  RenormalisedHOParameters coarse_param(param.M_lat,
                                        param.T_final,
                                        param.m0,
                                        param.mu2,
                                        param.perturbative);
  HarmonicOscillatorAction action(param.M_lat,
                                  param.T_final,
                                  param.m0,
                                  param.mu2);
  Sampler* sampler;
  if (param.hmc_sampling) {
    sampler = new HMCSampler(action,
                             param.T_hmc,
                             param.dt_hmc,
                             param.n_burnin_hmc);
  } else {
    sampler = &action;
  }
  HarmonicOscillatorAction coarse_action(param.M_lat/2,
                                         param.T_final,
                                         coarse_param.m0_coarse(),
                                         coarse_param.mu2_coarse());
  MonteCarloSingleLevel montecarlo_singlelevel(action,
                                               *sampler,
                                               qoi,
                                               param.n_samples,
                                               param.n_burnin);
  double exact_result = action.Xsquared_exact();
  double exact_result_continuum = action.Xsquared_exact_continuum();
  std::cout << std::endl;
  std::cout << " Exact result             <x^2> = " << exact_result << std::endl;
  std::cout << " Continuum limit [a -> 0] <x^2> = " << exact_result_continuum << std::endl;
  std::cout << std::endl;

  std::pair<double,double> result;
  result = montecarlo_singlelevel.evaluate();
  std::cout << " <x^2> = " << result.first << " +/- " << result.second << std::endl;
  std::cout << std::endl;
  Sampler* coarse_sampler;
  if (param.hmc_sampling) {
    coarse_sampler = new HMCSampler(coarse_action,
                                    param.T_hmc,
                                    param.dt_hmc,
                                    param.n_burnin_hmc);
  } else {
    coarse_sampler = &coarse_action;
  }
  MonteCarloTwoLevel montecarlo_twolevel(coarse_action,
                                         *coarse_sampler,
                                         action,
                                         qoi,
                                         param.n_samples,
                                         param.n_burnin,
                                         true);  
  result = montecarlo_twolevel.evaluate_difference();
  std::cout << " difference <x^2> " << std::endl;
  std::cout << " mean = " << result.first << " variance " << result.second << std::endl;
  std::cout << std::endl;
  std::cout << "=== Fine level sampler statistics === " << std::endl; 
  sampler->show_stats();
  std::cout << std::endl;
  std::cout << "=== Coarse level sampler statistics === " << std::endl; 
  coarse_sampler->show_stats();
    std::cout << std::endl;
  std::cout << "=== Two level sampler statistics === " << std::endl; 
  montecarlo_twolevel.get_twolevelsampler().show_stats();
  std::cout << std::endl;
  if (param.hmc_sampling) {
    delete sampler;
    delete coarse_sampler;
  }
}
