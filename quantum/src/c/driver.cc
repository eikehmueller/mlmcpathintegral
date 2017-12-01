#include <iostream>
#include <utility>
#include "harmonicoscillatoraction.hh"
#include "quantityofinterest.hh"
#include "montecarlo.hh"
#include "parameters.hh"

/** @file driver.cc
 * @brief Main program
 */

int main(int argc, char* argv[]) {
  std::cout << "*** Path integral multilevel MCMC ***" << std::endl;  
  Parameters param("parameters.in");
  param.show();
  QoIXsquared qoi(param.M_lat);
  HarmonicOscillatorAction action(param.M_lat,
                                  param.T_final,
                                  param.m0,
                                  param.mu2);
  HarmonicOscillatorAction coarse_action(param.M_lat/2,
                                         param.T_final,
                                         param.m0,
                                         param.mu2);
  MonteCarloSingleLevel montecarlo_singlelevel(action,
                                               action,
                                               qoi,
                                               param.n_samples,
                                               param.n_burnin);
  std::pair<double,double> result;
  double exact_result = action.Xsquared_exact();
  std::cout << std::endl;
  std::cout << " Exact result  <x^2> = " << exact_result << std::endl;
  std::cout << std::endl;
  result = montecarlo_singlelevel.evaluate();
  std::cout << " <x^2> = " << result.first << " +/- " << result.second << std::endl;
  std::cout << std::endl;
  
  MonteCarloTwoLevel montecarlo_twolevel(coarse_action,
                                         coarse_action,
                                         action,
                                         qoi,
                                         param.n_samples,
                                         param.n_burnin,
                                         true);  
  result = montecarlo_twolevel.evaluate_difference();
  std::cout << " difference <x^2> " << std::endl;
  std::cout << " mean = " << result.first << " variance " << result.second << std::endl;

  montecarlo_twolevel.show_stats();
}
