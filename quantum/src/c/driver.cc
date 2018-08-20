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
  std::cout << "   E[QoI] = " << result.first << " +/- " << error << std::endl;
  std::cout << " Var[QoI] = " << result.second << std::endl;
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
  result = montecarlo_twolevel.evaluate_difference();
  std::cout << " difference <QoI> " << std::endl;
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
