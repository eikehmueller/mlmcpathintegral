#include <iostream>
#include <utility>
#include "harmonicoscillatoraction.hh"
#include "quantityofinterest.hh"
#include "montecarlo.hh"

/** Main program
 */

int main(int argc, char* argv[]) {
  std::cout << "*** Path integral multilevel MCMC ***" << std::endl;  
  const unsigned int M_lat=16;
  const double T_final=1.0;
  const double m0=1.0;
  const double mu2=1.0;
  const unsigned int n_burnin=10000;
  const unsigned int n_samples=1000000;
  std::cout << " M_lat      = " << M_lat << std::endl;
  std::cout << " T_final    = " << T_final << std::endl;
  std::cout << " a_lat      = " << (T_final/M_lat) << std::endl;
  std::cout << " m_0        = " << m0 << std::endl;
  std::cout << " mu^2       = " << mu2 << std::endl;
  std::cout << " N_burnin   = " << n_burnin << std::endl;
  std::cout << " N_samples  = " << n_samples << std::endl;
  QoIXsquared qoi(M_lat);
  HarmonicOscillatorAction action(M_lat,T_final,m0,mu2);
  HarmonicOscillatorAction coarse_action(M_lat/2,T_final,m0,mu2);
  MonteCarloSingleLevel montecarlo_singlelevel(action,qoi,n_samples,n_burnin);
  std::pair<double,double> result;
  result = montecarlo_singlelevel.evaluate();
  std::cout << " <x^2> = " << result.first << " +/- " << result.second << std::endl;
  std::cout << std::endl;
  
  MonteCarloTwoLevel montecarlo_twolevel(coarse_action,
                                         coarse_action,
                                         action,
                                         qoi,
                                         n_samples,
                                         n_burnin);  
  result = montecarlo_twolevel.evaluate_difference();
  std::cout << " difference <x^2> " << std::endl;
  std::cout << " mean = " << result.first << " variance " << result.second << std::endl;
}
