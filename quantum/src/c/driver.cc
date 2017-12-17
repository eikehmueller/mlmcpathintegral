#include <iostream>
#include <utility>
#include <memory>
#include "harmonicoscillatoraction.hh"
#include "quarticoscillatoraction.hh"
#include "quantityofinterest.hh"
#include "montecarlo.hh"
#include "parameters.hh"
#include "renormalisation.hh"
#include "twolevelmetropolissampler.hh"
#include "hmcsampler.hh"
#include "config.h"

/** @file driver.cc
 * @brief File with main program
 *
 * @mainpage
 * Several classes for implementating Multilevel MCMC for the path-integral
 * formulation of quantum mechanics.
 */

/** Main program */
int main(int argc, char* argv[]) {
  std::cout << "+-----------------------------------+!" << std::endl;  
  std::cout << "!   Path integral multilevel MCMC   !" << std::endl;
  std::cout << "+-----------------------------------+!" << std::endl;  
  std::cout << std::endl;
  Parameters param("parameters.in");
  param.show();
  QoIXsquared qoi(param.M_lat);
  RenormalisedHOParameters coarse_param(param.M_lat,
                                        param.T_final,
                                        param.m0,
                                        param.mu2,
                                        param.perturbative);
  std::cout << std::endl;
  /* *** HARMONIC OSCILLATOR *** */
#ifdef ACTION_HARMONIC_OSCILLATOR
  HarmonicOscillatorAction action(param.M_lat,
                                  param.T_final,
                                  param.m0,
                                  param.mu2);
  HarmonicOscillatorAction coarse_action(param.M_lat/2,
                                         param.T_final,
                                         coarse_param.m0_coarse(),
                                         coarse_param.mu2_coarse());
  std::cout << "Action = harmonic oscillator" << std::endl;
#endif // ACTION_HARMONIC_OSCILLATOR
    /* *** QUARTIC OSCILLATOR *** */
#ifdef ACTION_QUARTIC_OSCILLATOR
  QuarticOscillatorAction action(param.M_lat,
                                 param.T_final,
                                 param.m0,
                                 param.mu2,
                                 param.lambda);
  QuarticOscillatorAction coarse_action(param.M_lat/2,
                                        param.T_final,
                                        param.m0,
                                        param.mu2,
                                        param.lambda);
  std::cout << "action = quartic oscillator" << std::endl;
#endif // QUARTIC_HARMONIC_OSCILLATOR
  std::cout << std::endl;
  
  Sampler* sampler;
  if (param.hmc_sampling) {
    sampler = new HMCSampler(action,
                             param.T_hmc,
                             param.dt_hmc,
                             param.n_burnin_hmc);
  } else {
#ifdef ACTION_HARMONIC_OSCILLATOR
    sampler = &action;
#else
    std::cout << " ERROR: can only sample directly from harmonic oscillator action." << std::endl;
    exit(-1);
#endif // ACTION_HARMONIC_OSCILLATOR
  }
  MonteCarloSingleLevel montecarlo_singlelevel(action,
                                               *sampler,
                                               qoi,
                                               param.n_samples,
                                               param.n_burnin);
#ifdef ACTION_HARMONIC_OSCILLATOR
  double exact_result = action.Xsquared_exact();
  double exact_result_continuum = action.Xsquared_exact_continuum();
  std::cout << std::endl;
  std::cout << " Exact result             <x^2> = " << exact_result << std::endl;
  std::cout << " Continuum limit [a -> 0] <x^2> = " << exact_result_continuum << std::endl;
  std::cout << std::endl;
#endif // ACTION_HARMONIC_OSCILLATOR
  std::pair<double,double> result;
  result = montecarlo_singlelevel.evaluate();
  double error = sqrt(result.second/(1.*param.n_samples));
  std::cout << "   E[x^2] = " << result.first << " +/- " << error << std::endl;
  std::cout << " Var[x^2] = " << result.second << std::endl;
  std::cout << std::endl;
  Sampler* coarse_sampler;
  if (param.hmc_sampling) {
    coarse_sampler = new HMCSampler(coarse_action,
                                    param.T_hmc,
                                    param.dt_hmc,
                                    param.n_burnin_hmc);
  } else {
#ifdef ACTION_HARMONIC_OSCILLATOR
    coarse_sampler = &coarse_action;
#else
    std::cout << " ERROR: can only sample directly from harmonic oscillator action." << std::endl;
    exit(-1);
#endif // ACTION_HARMONIC_OSCILLATOR
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

  // Tidy up
  if (param.hmc_sampling) {
    delete sampler;
    delete coarse_sampler;
  }
}
