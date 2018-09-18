#include <iostream>
#include <iomanip>
#include <utility>
#include <memory>
#include "harmonicoscillatoraction.hh"
#include "quarticoscillatoraction.hh"
#include "doublewellaction.hh"
#include "rotoraction.hh"
#include "quantityofinterest.hh"
#include "montecarlosinglelevel.hh"
#include "montecarlotwolevel.hh"
#include "montecarlomultilevel.hh"
#include "parameters.hh"
#include "renormalisation.hh"
#include "twolevelmetropolisstep.hh"
#include "conditionedfineaction.hh"
#include "mcmcstep.hh"
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
  LatticeParameters param_lattice;
  HarmonicOscillatorParameters param_ho;
  QuarticOscillatorParameters param_qo;
  DoubleWellParameters param_dw;
  RotorParameters param_rotor;
  if (param_general.readFile(filename)) return 1;
  std::cout << param_general << std::endl;
  if (param_lattice.readFile(filename)) return 1;
  std::cout << param_lattice << std::endl;
  switch (param_general.action()) {
  case (ActionHarmonicOscillator): {
    if (param_ho.readFile(filename)) return 1;
    std::cout << param_ho << std::endl;
    break;
  }
  case (ActionQuarticOscillator): {
    if (param_qo.readFile(filename)) return 1;
    std::cout << param_qo << std::endl;
    break;
  }
  case (ActionDoubleWell): {
    if (param_dw.readFile(filename)) return 1;
    std::cout << param_dw << std::endl;
    break;
  }
  case (ActionRotor): {
    if (param_rotor.readFile(filename)) return 1;
    std::cout << param_rotor << std::endl;
    break;
  }
  }

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

  MultiLevelMCParameters param_multilevelmc;
  if (param_multilevelmc.readFile(filename)) return 1;
  std::cout << param_multilevelmc << std::endl;

  /* ====== Select quantity of interest ====== */
  std::shared_ptr<QoI> qoi;
  std::cout << std::endl;
  if ( (param_general.action() == ActionHarmonicOscillator) or
       (param_general.action() == ActionQuarticOscillator) or
       (param_general.action() == ActionDoubleWell) ) {
    qoi=std::make_shared<QoIXsquared>();
    std::cout << "QoI = X^2 " << std::endl;
  }
  if ( (param_general.action() == ActionRotor) ) {
    qoi=std::make_shared<QoISusceptibility>();
    std::cout << "QoI = Susceptibility Q[X]^2/T " << std::endl;
  }
  std::cout << std::endl;
  
  /* ====== Select action ====== */
  std::shared_ptr<Action> action;
  switch (param_general.action()) {
  case (ActionHarmonicOscillator): {
    action =
      std::make_shared<HarmonicOscillatorAction>(param_lattice.M_lat(),
                                                 param_lattice.T_final(),
                                                 param_ho.renormalisation(),
                                                 param_ho.m0(),
                                                 param_ho.mu2());
    break;
  }
  case (ActionQuarticOscillator): {
    action =
      std::make_shared<QuarticOscillatorAction>(param_lattice.M_lat(),
                                                param_lattice.T_final(),
                                                RenormalisationNone,
                                                param_qo.m0(),
                                                param_qo.mu2(),
                                                param_qo.lambda());    
    break;
  }
  case (ActionDoubleWell): {
    action = 
      std::make_shared<DoubleWellAction>(param_lattice.M_lat(),
                                         param_lattice.T_final(),
                                         RenormalisationNone,
                                         param_dw.m0(),
                                         param_dw.mu2(),
                                         param_dw.lambda(),
                                         param_dw.sigma());

    break;
  }
  case (ActionRotor): {
    action = 
      std::make_shared<RotorAction>(param_lattice.M_lat(),
                                    param_lattice.T_final(),
                                    RenormalisationNone,
                                    param_rotor.m0());

    break;
  }
  }

  /* **************************************** * 
   * Single level method                      *
   * **************************************** */  

  if (param_general.do_singlelevelmc()) {
    std::cout << "+--------------------------------+" << std::endl;
    std::cout << "! Single level MC                !" << std::endl;
    std::cout << "+--------------------------------+" << std::endl;
    std::cout << std::endl;

    /* ====== Construct single level MC ====== */
    MonteCarloSingleLevel montecarlo_singlelevel(action,
                                                 qoi,
                                                 param_general,
                                                 param_hmc,
                                                 param_cluster,
                                                 param_singlelevelmc);

    /* ====== Print out exact result for harmonic oscillator */
    if (param_general.action() == ActionHarmonicOscillator) {
      std::shared_ptr<HarmonicOscillatorAction> ho_action =
        std::dynamic_pointer_cast<HarmonicOscillatorAction>(action);
      double exact_result = ho_action->Xsquared_exact();
      double exact_result_continuum = ho_action->Xsquared_exact_continuum();
      std::cout << std::endl;
      std::cout << std::setprecision(6) << std::fixed;
      std::cout << " Exact result             <x^2> = " << exact_result << std::endl;
      std::cout << " Continuum limit [a -> 0] <x^2> = " << exact_result_continuum << std::endl;
      std::cout << std::endl;
    }
  
    Statistics stats("QoI",10);
    montecarlo_singlelevel.evaluate(stats);
    std::cout << stats << std::endl;
    std::cout << "=== Sampler statistics === " << std::endl; 
    montecarlo_singlelevel.get_sampler()->show_stats();
    std::cout << std::endl;

  }

  /* **************************************** * 
   * Two level method                         *
   * **************************************** */
  if (param_general.do_twolevelmc()) {
    std::cout << "+--------------------------------+" << std::endl;
    std::cout << "! Two level MC                   !" << std::endl;
    std::cout << "+--------------------------------+" << std::endl;
    std::cout << std::endl;
    
    MonteCarloTwoLevel montecarlo_twolevel(action,
                                           qoi,
                                           param_general,
                                           param_hmc,
                                           param_cluster,
                                           param_twolevelmc);
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
    std::cout << "=== Coarse level sampler statistics === " << std::endl; 
    montecarlo_twolevel.get_coarsesampler()->show_stats();
    std::cout << std::endl;
    std::cout << "=== Two level sampler statistics === " << std::endl; 
    montecarlo_twolevel.get_twolevelstep()->show_stats();
    std::cout << std::endl;
  }
  
  /* **************************************** * 
   * Multilevel method                         *
   * **************************************** */
  if (param_general.do_multilevelmc()) {
    std::cout << "+--------------------------------+" << std::endl;
    std::cout << "! Multilevel MC                  !" << std::endl;
    std::cout << "+--------------------------------+" << std::endl;
    std::cout << std::endl;
    
    MonteCarloMultiLevel montecarlo_multilevel(action,
                                               qoi,
                                               param_general,
                                               param_hmc,
                                               param_cluster,
                                               param_multilevelmc);
    montecarlo_multilevel.evaluate();
  }
}
