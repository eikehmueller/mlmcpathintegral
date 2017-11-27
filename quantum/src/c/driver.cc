#include <iostream>
#include <utility>
#include "harmonicoscillatoraction.hh"
#include "quantityofinterest.hh"
#include "montecarlo.hh"

/** Main program
 */

int main(int argc, char* argv[]) {
  std::cout << "*** Path integral multilevel MCMC ***" << std::endl;  
  const unsigned int M_lat=8;
  const double T_final=1.0;
  const double m0=1.0;
  const double mu2=1.0;
  const unsigned int N_burnin=10000;
  const unsigned int N_samples=100000;
  std::cout << " M_lat      = " << M_lat << std::endl;
  std::cout << " T_final    = " << T_final << std::endl;
  std::cout << " a_lat      = " << (T_final/M_lat) << std::endl;
  std::cout << " m_0        = " << m0 << std::endl;
  std::cout << " mu^2       = " << mu2 << std::endl;
  std::cout << " N_burnin   = " << N_burnin << std::endl;
  std::cout << " N_samples  = " << N_samples << std::endl;
  QoIXsquared qoi(M_lat);
  HarmonicOscillatorAction action(M_lat,T_final,m0,mu2);
  MonteCarloSingleLevel montecarlo_singlelevel(action,qoi,N_samples,N_burnin);
  std::pair<double,double> result = montecarlo_singlelevel.evaluate();
  std::cout << " <x^2> = " << result.first << " +/- " << result.second << std::endl;
}
