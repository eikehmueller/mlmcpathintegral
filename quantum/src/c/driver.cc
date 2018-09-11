#include <iostream>
#include <iomanip>
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
  std::cout << "++===================================++" << std::endl;
  std::cout << "!!   Path integral multilevel MCMC   !!" << std::endl;
  std::cout << "++===================================++" << std::endl;
  std::cout << std::endl;
  if (argc != 2) {
    std::cout << "Usage: " << argv[0] << " PARAMETERFILE" << std::endl;
    std::cout << std::endl;
    return 0;
  }
  std::string filename = argv[1];
  std::cout << " Reading parameter from file \'" << filename << "\'" << std::endl;
    std::cout << std::endl;
  
  /* ====== Read parameters ====== */
  GeneralParameters param_general;
  if (param_general.readFile(filename)) return 1;
  std::cout << param_general << std::endl;
  LatticeParameters param_lattice;
  if (param_lattice.readFile(filename)) return 1;
  std::cout << param_lattice << std::endl;
  
#ifdef ACTION_HARMONIC_OSCILLATOR
  HarmonicOscillatorParameters param_ho;
  if (param_ho.readFile(filename)) return 1;
  std::cout << param_ho << std::endl;
#endif // ACTION_HARMONIC_OSCILLATOR

#ifdef ACTION_QUARTIC_OSCILLATOR
  QuarticOscillatorParameters param_qo;
  if (param_qo.readFile(filename)) return 1;
  std::cout << param_qo << std::endl;
#endif // ACTION_QUARTIC_OSCILLATOR

#ifdef ACTION_DOUBLE_WELL
  DoubleWellParameters param_dw;
  if (param_dw.readFile(filename)) return 1;
  std::cout << param_dw << std::endl;
#endif // ACTION_DOUBLE_WELL

#ifdef ACTION_ROTOR
  RotorParameters param_rotor;
  if (param_rotor.readFile(filename)) return 1;
  std::cout << param_rotor << std::endl;
#endif // ACTION_ROTOR

  HMCParameters param_hmc;
  if (param_hmc.readFile(filename)) return 1;
  std::cout << param_hmc << std::endl;

  ClusterParameters param_cluster;
  if (param_cluster.readFile(filename)) return 1;
  std::cout << param_cluster << std::endl;

  SingleLevelMCParameters param_singlelevelmc;
  if (param_singlelevelmc.readFile(filename)) return 1;
  std::cout << param_singlelevelmc << std::endl;

  TwoLevelMCParameters param_twolevelmc;
  if (param_twolevelmc.readFile(filename)) return 1;
  std::cout << param_twolevelmc << std::endl;
  
  /* ====== Select quantity of interest ====== */
  std::shared_ptr<QoI> qoi;
#if defined(ACTION_HARMONIC_OSCILLATOR) || defined(ACTION_QUARTIC_OSCILLATOR) || defined(ACTION_DOUBLE_WELL)
  qoi=std::make_shared<QoIXsquared>(param_lattice.M_lat());
  std::cout << std::endl;
  std::cout << "QoI = X^2 " << std::endl;
  std::cout << std::endl;
#endif // ACTION_HARMONIC_OSCILLATOR || ACTION_QUARTIC_OSCILLATOR || ACTION_DOUBLEWELL 
#ifdef ACTION_ROTOR
  qoi=std::make_shared<QoISusceptibility>(param_lattice.M_lat());
  std::cout << std::endl;
  std::cout << "QoI = Susceptibility Q[X]^2/T " << std::endl;
  std::cout << std::endl;
#endif

  /* ====== Select action ====== */
  /* *** HARMONIC OSCILLATOR *** */
#ifdef ACTION_HARMONIC_OSCILLATOR
  std::cout << "Action = harmonic oscillator" << std::endl;
  std::shared_ptr<Action> action =
    std::make_shared<HarmonicOscillatorAction>(param_lattice.M_lat(),
                                               param_lattice.T_final(),
                                               param_ho.renormalisation(),
                                               param_ho.m0(),
                                               param_ho.mu2());
#endif // ACTION_HARMONIC_OSCILLATOR
  
  /* *** QUARTIC OSCILLATOR *** */
#ifdef ACTION_QUARTIC_OSCILLATOR
  std::cout << "Action = quartic oscillator" << std::endl;
  std::shared_ptr<Action> action =
    std::make_shared<QuarticOscillatorAction>(param_lattice.M_lat(),
                                              param_lattice.T_final(),
                                              RenormalisationNone,
                                              param_qo.m0(),
                                              param_qo.mu2(),
                                              param_qo.lambda());
#endif // ACTION_QUARTIC_OSCILLATOR
  
  /* *** DOUBLE WELL *** */
#ifdef ACTION_DOUBLE_WELL
  std::cout << "Action = double well" << std::endl;
  std::shared_ptr<Action> action = 
    std::make_shared<DoubleWellAction>(param_lattice.M_lat(),
                                       param_lattice.T_final(),
                                       RenormalisationNone,
                                       param_dw.m0(),
                                       param_dw.mu2(),
                                       param_dw.lambda(),
                                       param_dw.sigma());
#endif // ACTION_DOUBLE_WELL

  /* *** ROTOR *** */
#ifdef ACTION_ROTOR
  std::cout << "Action = rotor" << std::endl;
  std::shared_ptr<Action> action = 
    std::make_shared<RotorAction>(param_lattice.M_lat(),
                                  param_lattice.T_final(),
                                  RenormalisationNone,
                                  param_rotor.m0());
#endif // ACTION_ROTOR
  std::cout << std::endl;

  /* **************************************** * 
   * Single level method                      *
   * **************************************** */  
  /* ====== Select sampler for single level method ====== */
  std::shared_ptr<Sampler> sampler;
  if (param_singlelevelmc.sampler() == SamplerHMC) {
    sampler = std::make_shared<HMCSampler>(action,
                                           param_hmc.T(),
                                           param_hmc.dt(),
                                           param_hmc.n_burnin());
  } else if (param_singlelevelmc.sampler() == SamplerCluster) {
#ifdef ACTION_ROTOR
    sampler = std::make_shared<ClusterSampler>(std::dynamic_pointer_cast<ClusterAction>(action),
                                               param_cluster.n_burnin(),
                                               param_cluster.n_updates());
#else
    std::cout << " ERROR: can only use cluster sampler for QM rotor action." << std::endl;
    exit(-1);
#endif
  } else if (param_singlelevelmc.sampler() == SamplerExact) {
#ifdef ACTION_HARMONIC_OSCILLATOR
    sampler = std::dynamic_pointer_cast<Sampler>(action);
#else
    std::cout << " ERROR: can only sample exactly from harmonic oscillator action." << std::endl;
    exit(-1);
#endif // ACTION_HARMONIC_OSCILLATOR
  } else {
    std::cout << " ERROR: Unknown sampler." << std::endl;
    exit(-1);
  }
  if (param_general.do_singlelevelmc()) {
    std::cout << "+--------------------------------+" << std::endl;
    std::cout << "! Single level MC                !" << std::endl;
    std::cout << "+--------------------------------+" << std::endl;
    std::cout << std::endl;

    /* ====== Construct single level MC ====== */
    MonteCarloSingleLevel montecarlo_singlelevel(action,
                                                 sampler,
                                                 qoi,
                                                 param_singlelevelmc.n_samples(),
                                                 param_singlelevelmc.n_burnin());

    /* ====== Print out exact result for harmonic oscillator */
#ifdef ACTION_HARMONIC_OSCILLATOR
    std::shared_ptr<HarmonicOscillatorAction> ho_action =
      std::dynamic_pointer_cast<HarmonicOscillatorAction>(action);
    double exact_result = ho_action->Xsquared_exact();
    double exact_result_continuum = ho_action->Xsquared_exact_continuum();
    std::cout << std::endl;
    std::cout << std::setprecision(6) << std::fixed;
    std::cout << " Exact result             <x^2> = " << exact_result << std::endl;
    std::cout << " Continuum limit [a -> 0] <x^2> = " << exact_result_continuum << std::endl;
    std::cout << std::endl;
#endif // ACTION_HARMONIC_OSCILLATOR
  
    Statistics stats("QoI",10);
    montecarlo_singlelevel.evaluate(stats);
    std::cout << stats << std::endl;
  }
  /* **************************************** * 
   * Two level method                         *
   * **************************************** */
  /* ====== Select coarse action ====== */
  std::shared_ptr<Action> coarse_action;
  coarse_action = action->coarse_action();
  
  std::shared_ptr<Sampler> coarse_sampler;
  if (param_twolevelmc.coarsesampler() == SamplerHMC) {
    coarse_sampler = std::make_shared<HMCSampler>(coarse_action,
                                                  param_hmc.T(),
                                                  param_hmc.dt(),
                                                  param_hmc.n_burnin());
  } else if (param_twolevelmc.coarsesampler() == SamplerCluster) {
#ifdef ACTION_ROTOR
    coarse_sampler = std::make_shared<ClusterSampler>(std::dynamic_pointer_cast<ClusterAction>(coarse_action),
                                                      param_cluster.n_burnin(),
                                                      param_cluster.n_updates());
#else
    std::cout << " ERROR: can only use cluster sampler for QM rotor action." << std::endl;
    exit(-1);
#endif // ACTION_ROTOR
  } else if (param_twolevelmc.coarsesampler() == SamplerExact) {
#ifdef ACTION_HARMONIC_OSCILLATOR
    coarse_sampler = std::dynamic_pointer_cast<Sampler>(coarse_action);
    
#else
    std::cout << " ERROR: can only sample exactly from harmonic oscillator action." << std::endl;
    exit(-1);
#endif // ACTION_HARMONIC_OSCILLATOR
  } else {
    std::cout << " ERROR: Unknown sampler." << std::endl;
    exit(-1);
  }

  std::shared_ptr<ConditionedFineAction> conditioned_fine_action;
#ifdef ACTION_ROTOR
  conditioned_fine_action =
    std::make_shared<RotorConditionedFineAction>(std::dynamic_pointer_cast<RotorAction>(action));
#else
  conditioned_fine_action =
    std::make_shared<GaussianConditionedFineAction>(action);
                                                                            
#endif // ACTION_ROTOR
  if (param_general.do_twolevelmc()) {
    std::cout << "+--------------------------------+" << std::endl;
    std::cout << "! Two level MC                   !" << std::endl;
    std::cout << "+--------------------------------+" << std::endl;
    std::cout << std::endl;
    MonteCarloTwoLevel montecarlo_twolevel(coarse_action,
                                           coarse_sampler,
                                           action,
                                           conditioned_fine_action,
                                           qoi,
                                           param_twolevelmc.n_samples(),
                                           param_twolevelmc.n_burnin(),
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
}
