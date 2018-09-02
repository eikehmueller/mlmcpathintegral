#include <iostream>
#include <utility>
#include <memory>
#include "harmonicoscillatoraction.hh"
#include "quarticoscillatoraction.hh"
#include "doublewellaction.hh"
#include "rotoraction.hh"
#include "quantityofinterest.hh"
#include "montecarlo.hh"
#include "parameters.hh"
#include "renormalisation.hh"
#include "twolevelmetropolissampler.hh"
#include "conditionedfineaction.hh"
#include "hmcsampler.hh"
#include "clustersampler.hh"
#include "config.h"
#include "statistics.hh"

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
#if defined(ACTION_HARMONIC_OSCILLATOR) || defined(ACTION_QUARTIC_OSCILLATOR) || defined(ACTION_DOUBLE_WELL)
  QoIXsquared qoi(param.M_lat);
  std::cout << std::endl;
  std::cout << "QoI = X^2 " << std::endl;
  std::cout << std::endl;
#endif // ACTION_HARMONIC_OSCILLATOR || ACTION_QUARTIC_OSCILLATOR || ACTION_DOUBLEWELL 
#ifdef ACTION_ROTOR
  QoISusceptibility qoi(param.M_lat);
  std::cout << std::endl;
  std::cout << "QoI = Susceptibility Q[X]^2/T " << std::endl;
  std::cout << std::endl;
#endif
  RenormalisedHOParameters coarse_param(param.M_lat,
                                        param.T_final,
                                        param.m0,
                                        param.mu2,
                                        param.renormalisation);
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
  std::cout << "Action = quartic oscillator" << std::endl;
#endif // ACTION_QUARTIC_OSCILLATOR
  /* *** DOUBLE WELL *** */
#ifdef ACTION_DOUBLE_WELL
  DoubleWellAction action(param.M_lat,
                          param.T_final,
                          param.m0,
                          param.mu2,
                          param.lambda,
                          param.sigma);
  DoubleWellAction coarse_action(param.M_lat/2,
                                 param.T_final,
                                 param.m0,
                                 param.mu2,
                                 param.lambda,
                                 param.sigma);
  std::cout << "Action = double well" << std::endl;
#endif // ACTION_DOUBLE_WELL
  /* *** ROTOR *** */
#ifdef ACTION_ROTOR
  RotorAction action(param.M_lat,
                     param.T_final,
                     param.m0);
  RotorAction coarse_action(param.M_lat/2,
                            param.T_final,
                            param.m0);
  std::cout << "Action = rotor" << std::endl;
#endif // ACTION_ROTOR
  std::cout << std::endl;
  
  std::shared_ptr<Sampler> sampler;
  if (param.sampler == SamplerHMC) {
    sampler = std::make_shared<HMCSampler>(action,
                                           param.T_hmc,
                                           param.dt_hmc,
                                           param.n_burnin_sampler);
  } else if (param.sampler == SamplerCluster) {
#ifdef ACTION_ROTOR
    sampler = std::make_shared<ClusterSampler>(action,
                                               param.n_burnin_sampler);
#else
    std::cout << " ERROR: can only use cluster sampler for QM rotor action." << std::endl;
#endif
  } else if (param.sampler == SamplerExact) {
#ifdef ACTION_HARMONIC_OSCILLATOR
    sampler = &action;
#else
    std::cout << " ERROR: can only sample exactly from harmonic oscillator action." << std::endl;
    exit(-1);
#endif // ACTION_HARMONIC_OSCILLATOR
  } else {
    std::cout << " ERROR: Unknown sampler." << std::endl;
    exit(-1);
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
  Statistics stats("QoI",10);
  montecarlo_singlelevel.evaluate(stats);
  std::cout << stats << std::endl;
  std::shared_ptr<Sampler> coarse_sampler;
  if (param.sampler == SamplerHMC) {
    coarse_sampler = std::make_shared<HMCSampler>(coarse_action,
                                                  param.T_hmc,
                                                  param.dt_hmc,
                                                  param.n_burnin_sampler);
  } else if (param.sampler == SamplerCluster) {
#ifdef ACTION_ROTOR
    coarse_sampler = std::make_shared<ClusterSampler>(coarse_action,
                                                      param.n_burnin_sampler);
#else
    std::cout << " ERROR: can only use cluster sampler for QM rotor action." << std::endl;
    exit(-1);
#endif // ACTION_ROTOR
  } else if (param.sampler == SamplerExact) {
#ifdef ACTION_HARMONIC_OSCILLATOR
    coarse_sampler = &coarse_action;
#else
    std::cout << " ERROR: can only sample exactly from harmonic oscillator action." << std::endl;
    exit(-1);
#endif // ACTION_HARMONIC_OSCILLATOR
  } else {
    std::cout << " ERROR: Unknown sampler." << std::endl;
    exit(-1);
  }
#ifdef ACTION_ROTOR
  RotorConditionedFineAction conditioned_fine_action(action);
#else
  GaussianConditionedFineAction conditioned_fine_action(action);
#endif // ACTION_ROTOR
  // Two level method
  MonteCarloTwoLevel montecarlo_twolevel(coarse_action,
                                         *coarse_sampler,
                                         action,
                                         conditioned_fine_action,
                                         qoi,
                                         param.n_samples,
                                         param.n_burnin,
                                         true);
  Statistics stats_fine("QoI[fine]",10);
  Statistics stats_coarse("QoI[coarse]",10);
  Statistics stats_diff("delta QoI",10);
  montecarlo_twolevel.evaluate_difference(stats_fine,
                                          stats_coarse,
                                          stats_diff);
  std::cout << stats_fine << std::endl;
  std::cout << stats_coarse << std::endl;
  std::cout << stats_diff << std::endl;
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
}
