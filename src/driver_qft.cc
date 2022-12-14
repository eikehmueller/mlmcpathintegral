#include <iomanip>
#include <iostream>
#include <memory>
#include <utility>

#include "action/qft/gffaction.hh"
#include "action/qft/gffconditionedfineaction.hh"
#include "action/qft/nonlinearsigmaaction.hh"
#include "action/qft/nonlinearsigmaconditionedfineaction.hh"
#include "action/qft/quenchedschwingeraction.hh"
#include "action/qft/quenchedschwingerconditionedfineaction.hh"
#include "action/renormalisation.hh"
#include "common/parameters.hh"
#include "common/statistics.hh"
#include "config.h"
#include "lattice/lattice2d.hh"
#include "montecarlo/montecarlomultilevel.hh"
#include "montecarlo/montecarlosinglelevel.hh"
#include "montecarlo/montecarlotwolevel.hh"
#include "mpi/mpi_wrapper.hh"
#include "qoi/qft/qoi2dmagneticsusceptibility.hh"
#include "qoi/qft/qoi2dphisquared.hh"
#include "qoi/qft/qoi2dsusceptibility.hh"
#include "qoi/qft/qoiavgplaquette.hh"
#include "sampler/clustersampler.hh"
#include "sampler/hierarchicalsampler.hh"
#include "sampler/hmcsampler.hh"
#include "sampler/multilevelsampler.hh"
#include "sampler/overrelaxedheatbathsampler.hh"
#include "sampler/quenchedschwingerclustersampler.hh"

/** @file driver_schwinger.cc
 * @brief File with main program for 2D Schwinger model
 *
 * @mainpage
 * Several classes for implementating Multilevel MCMC for the path-integral
 * formulation of the lattice Schwinger model in 2D.
 */

/** Helper function to construct suitable sampler factory for given samplerid */
std::shared_ptr<SamplerFactory> construct_sampler_factory(
    const int samplerid, const std::shared_ptr<QoIFactory> qoi_factory,
    const std::shared_ptr<SamplerFactory> coarse_sampler_factory,
    const std::shared_ptr<ConditionedFineActionFactory>
        conditioned_fine_action_factory,
    const GeneralParameters param_general, const QFTParameters param_qft,
    const HMCParameters param_hmc, const ClusterParameters param_cluster,
    const OverrelaxedHeatBathParameters param_heatbath,
    const HierarchicalParameters param_hierarchical,
    const StatisticsParameters param_stats) {
  std::shared_ptr<SamplerFactory> sampler_factory;
  if (samplerid == SamplerHMC) {
    /* --- CASE 1: HMC sampler ---- */
    sampler_factory = std::make_shared<HMCSamplerFactory>(param_hmc);
  } else if (samplerid == SamplerOverrelaxedHeatBath) {
    /* --- CASE 2: heat bath sampler ---- */
    sampler_factory =
        std::make_shared<OverrelaxedHeatBathSamplerFactory>(param_heatbath);
  } else if (samplerid == SamplerHierarchical) {
    /* --- CASE 3: Hierarchical sampler */
    sampler_factory = std::make_shared<HierarchicalSamplerFactory>(
        coarse_sampler_factory, conditioned_fine_action_factory,
        param_hierarchical);
  } else if (samplerid == SamplerMultilevel) {
    /* --- CASE 4: Multilevel sampler */
    sampler_factory = std::make_shared<MultilevelSamplerFactory>(
        qoi_factory, coarse_sampler_factory, conditioned_fine_action_factory,
        param_stats, param_hierarchical);
  } else if (samplerid == SamplerCluster) {
    /* --- CASE 5: Cluster sampler */
    if (param_qft.action() == ActionQuenchedSchwinger) {
      sampler_factory =
          std::make_shared<QuenchedSchwingerClusterSamplerFactory>(
              param_cluster);
    } else if (param_qft.action() == ActionNonlinearSigma) {
      sampler_factory = std::make_shared<ClusterSamplerFactory>(param_cluster);
    } else {
      mpi_parallel::cerr << " ERROR: cluster not supported for chosen action."
                         << std::endl;
      mpi_exit(EXIT_FAILURE);
    }
  } else if (samplerid == SamplerExact) {
    /* --- CASE 6: Exact sampler */
    if (param_qft.action() == ActionGFF) {
      sampler_factory = std::make_shared<GFFSamplerFactory>();
    } else {
      mpi_parallel::cerr
          << " ERROR: exact sampler not supported for chosen action."
          << std::endl;
      mpi_exit(EXIT_FAILURE);
    }
  } else {
    mpi_parallel::cerr << " ERROR: Unsupported sampler." << std::endl;
    mpi_exit(EXIT_FAILURE);
  }
  return sampler_factory;
}

/** Main program */
int main(int argc, char *argv[]) {
  mpi_init();
  Timer total_time("total");
  total_time.start();
  mpi_parallel::cout << "++===================================++" << std::endl;
  mpi_parallel::cout << "!!   Path integral multilevel MCMC   !!" << std::endl;
  mpi_parallel::cout << "!!   for the 2D Schwinger model      !!" << std::endl;
  mpi_parallel::cout << "++===================================++" << std::endl;
  mpi_parallel::cout << std::endl;
#ifdef USE_MPI
  mpi_parallel::cout << "MPI parallel version running on " << mpi_comm_size()
                     << " processes." << std::endl;
#else
  mpi_parallel::cout << "Sequential version." << std::endl;
#endif // USE_MPI
  mpi_parallel::cout << std::endl;
  mpi_parallel::cout << "Starting run at " << current_time() << std::endl;
  if (argc != 2) {
    mpi_parallel::cout << "Usage: " << argv[0] << " PARAMETERFILE" << std::endl;
    mpi_parallel::cout << std::endl;
    return 0;
  }
  std::string filename = argv[1];
  mpi_parallel::cout << " Reading parameter from file \'" << filename << "\'"
                     << std::endl;
  mpi_parallel::cout << std::endl;

  /* ====== Read parameters ====== */

  GeneralParameters param_general;
  if (param_general.readFile(filename))
    return 1;
  mpi_parallel::cout << param_general << std::endl;

  QFTParameters param_qft;
  if (param_qft.readFile(filename))
    return 1;
  mpi_parallel::cout << param_qft << std::endl;

  Lattice2DParameters param_lattice;
  if (param_lattice.readFile(filename))
    return 1;
  mpi_parallel::cout << param_lattice << std::endl;

  StatisticsParameters param_stats;
  if (param_stats.readFile(filename))
    return 1;
  mpi_parallel::cout << param_stats << std::endl;

  NonlinearSigmaParameters param_nonlinearsigma;
  SchwingerParameters param_schwinger;
  GFFParameters param_gff;
  switch (param_qft.action()) {
  case (ActionQuenchedSchwinger): {
    if (param_schwinger.readFile(filename))
      return 1;
    mpi_parallel::cout << param_schwinger << std::endl;
    break;
  }
  case (ActionNonlinearSigma): {
    if (param_nonlinearsigma.readFile(filename))
      return 1;
    mpi_parallel::cout << param_nonlinearsigma << std::endl;
    break;
  }
  case (ActionGFF): {
    if (param_gff.readFile(filename))
      return 1;
    mpi_parallel::cout << param_gff << std::endl;
    break;
  }
  }

  HMCParameters param_hmc;
  if (param_hmc.readFile(filename))
    return 1;
  mpi_parallel::cout << param_hmc << std::endl;

  ClusterParameters param_cluster;
  if (param_cluster.readFile(filename))
    return 1;
  mpi_parallel::cout << param_cluster << std::endl;

  OverrelaxedHeatBathParameters param_heatbath;
  if (param_heatbath.readFile(filename))
    return 1;
  mpi_parallel::cout << param_heatbath << std::endl;

  SingleLevelMCParameters param_singlelevelmc;
  if (param_singlelevelmc.readFile(filename))
    return 1;
  mpi_parallel::cout << param_singlelevelmc << std::endl;

  HierarchicalParameters param_hierarchical;
  if (param_hierarchical.readFile(filename))
    return 1;
  mpi_parallel::cout << param_hierarchical << std::endl;

  TwoLevelMCParameters param_twolevelmc;
  if (param_twolevelmc.readFile(filename))
    return 1;
  mpi_parallel::cout << param_twolevelmc << std::endl;

  MultiLevelMCParameters param_multilevelmc;
  if (param_multilevelmc.readFile(filename))
    return 1;
  mpi_parallel::cout << param_multilevelmc << std::endl;

#ifdef DEBUG_BUILD
  mpi_parallel::cout << FRED("CAUTION: built in debug mode.") << std::endl;
#endif // DEBUG_BUILD

#ifdef OPT_BUILD
  mpi_parallel::cout << FGREEN("Built in optimised mode.") << std::endl;
#endif // OPT_BUILD

#ifdef SAVE_PATHS
  mpi_parallel::cout << FRED("CAUTION: logging paths will impact performance!")
                     << std::endl;
#endif // SAVE_PATHS

#ifdef LOG_QOI
  mpi_parallel::cout << FRED("CAUTION: logging QoI will impact performance!")
                     << std::endl;
#endif // LOG_QOI

  /* ====== Lattice ====== */
  std::shared_ptr<Lattice2D> lattice;
  lattice = std::make_shared<Lattice2D>(param_lattice.Mt_lat(),
                                        param_lattice.Mx_lat(),
                                        param_lattice.coarsening_type());

  /* ====== Select quantity of interest and construct QoI factory ====== */
  std::shared_ptr<QoI> qoi;
  std::shared_ptr<QoIFactory> qoi_factory;
  mpi_parallel::cout << std::endl;
  if (param_qft.action() == ActionQuenchedSchwinger) {
    qoi = std::make_shared<QoI2DSusceptibility>(lattice);
    qoi_factory = std::make_shared<QoI2DSusceptibilityFactory>();
    mpi_parallel::cout << "QoI = Susceptibility Q[phi]^2 " << std::endl;
  }
  if ((param_qft.action() == ActionNonlinearSigma)) {
    qoi = std::make_shared<QoI2DMagneticSusceptibility>(lattice);
    qoi_factory = std::make_shared<QoI2DMagneticSusceptibilityFactory>();
    mpi_parallel::cout << "QoI = Average squared magnetisation 1/M*mu[phi]^2 "
                       << std::endl;
  }
  if ((param_qft.action() == ActionGFF)) {
    qoi = std::make_shared<QoI2DPhiSquared>(lattice);
    qoi_factory = std::make_shared<QoI2DPhiSquaredFactory>();
    mpi_parallel::cout << "QoI = Mean squared field 1/M*sum phi^2 "
                       << std::endl;
  }

  /* ====== Select action ====== */
  std::shared_ptr<Action> action;
  switch (param_qft.action()) {
  case (ActionQuenchedSchwinger): {
    action = std::make_shared<QuenchedSchwingerAction>(
        lattice, nullptr, param_schwinger.renormalisation(),
        param_schwinger.beta());
    break;
  }
  case (ActionNonlinearSigma): {
    action = std::make_shared<NonlinearSigmaAction>(
        lattice, nullptr, param_nonlinearsigma.renormalisation(),
        param_nonlinearsigma.beta());
    break;
  }
  case (ActionGFF): {
    action = std::make_shared<GFFAction>(lattice, nullptr, param_gff.mass());
    break;
  }
  }

  // numerical result, statistical error and analytical result
  double numerical_result;
  double statistical_error;
  double analytical_result;

  // do we compare to analytical result (only supported for some actions)
  bool do_analytical_comparison = false;
  if (param_qft.action() == ActionQuenchedSchwinger) {
    do_analytical_comparison = true;
    if (param_schwinger.beta() > 2000.0) {
      analytical_result = quenchedschwinger_chit_perturbative(
          param_schwinger.beta(), lattice->getNcells());
    } else {
      analytical_result = quenchedschwinger_chit_analytical(
          param_schwinger.beta(), lattice->getNcells());
    }
    double analytical_result_variance =
        quenchedschwinger_var_chit_continuum_analytical(param_schwinger.beta(),
                                                        lattice->getNcells());
    mpi_parallel::cout << std::endl;
    mpi_parallel::cout << std::setprecision(8) << std::fixed;
    mpi_parallel::cout << " Analytical results" << std::endl;
    mpi_parallel::cout << "      E[V*chi_t]              = "
                       << analytical_result;
    if (param_schwinger.beta() > 2000.0) {
      mpi_parallel::cout << " + O(beta^{-2}) = O("
                         << pow(param_schwinger.beta(), -2) << ")";
    }
    mpi_parallel::cout << std::endl;
    mpi_parallel::cout << "      lim_{a->0} Var[V*chi_t] = "
                       << analytical_result_variance << std::endl;
    mpi_parallel::cout << std::endl;
  } else if (param_qft.action() == ActionGFF) {
    do_analytical_comparison = true;
    analytical_result = gff_phi_squared_analytical(
        param_gff.mass(), lattice->getMt_lat(), lattice->getMx_lat());
    mpi_parallel::cout << std::endl;
    mpi_parallel::cout << std::setprecision(8) << std::fixed;
    mpi_parallel::cout << " Analytical result" << std::endl;
    mpi_parallel::cout << "      E[Q^2]              = " << analytical_result;
    mpi_parallel::cout << std::endl;
  }

  /* Construction conditioned fine action factory */
  std::shared_ptr<ConditionedFineActionFactory> conditioned_fine_action_factory;
  switch (param_qft.action()) {
  case (ActionQuenchedSchwinger): {
    conditioned_fine_action_factory =
        std::make_shared<QuenchedSchwingerConditionedFineActionFactory>();
    break;
  }
  case (ActionNonlinearSigma): {
    conditioned_fine_action_factory =
        std::make_shared<NonlinearSigmaConditionedFineActionFactory>();
    break;
  }
  case (ActionGFF): {
    conditioned_fine_action_factory =
        std::make_shared<GFFConditionedFineActionFactory>();
    break;
  }
  }

  /* Construct coarse level sampler factory, which might be used by the
   * hierarchical samplers */
  std::shared_ptr<SamplerFactory> coarse_sampler_factory;
  /* Note that here it does not make sense to use the hierarchical- or
   * multilevel-sampler, so we can pass null pointers for the QoI and coarse
   * level sampler factory
   */
  coarse_sampler_factory = construct_sampler_factory(
      param_hierarchical.coarsesampler(), nullptr, nullptr, nullptr,
      param_general, param_qft, param_hmc, param_cluster, param_heatbath,
      param_hierarchical, param_stats);

  /* **************************************** *
   * Single level method                      *
   * **************************************** */

  if (param_general.method() == MethodSingleLevel) {
    mpi_parallel::cout << "+--------------------------------+" << std::endl;
    mpi_parallel::cout << "! Single level MC                !" << std::endl;
    mpi_parallel::cout << "+--------------------------------+" << std::endl;
    mpi_parallel::cout << std::endl;

    std::shared_ptr<SamplerFactory> sampler_factory;
    sampler_factory = construct_sampler_factory(
        param_singlelevelmc.sampler(), qoi_factory, coarse_sampler_factory,
        conditioned_fine_action_factory, param_general, param_qft, param_hmc,
        param_cluster, param_heatbath, param_hierarchical, param_stats);

    /* ====== Construct single level MC ====== */
    MonteCarloSingleLevel montecarlo_singlelevel(
        action, qoi, sampler_factory, param_stats, param_singlelevelmc);

    montecarlo_singlelevel.evaluate();
    mpi_parallel::cout << std::endl;
    montecarlo_singlelevel.show_statistics();
    numerical_result = montecarlo_singlelevel.numerical_result();
    statistical_error = montecarlo_singlelevel.statistical_error();
    mpi_parallel::cout << "=== Sampler statistics === " << std::endl;
    montecarlo_singlelevel.get_sampler()->show_stats();
    mpi_parallel::cout << std::endl;
  }

  /* **************************************** *
   * Two level method                         *
   * **************************************** */
  if (param_general.method() == MethodTwoLevel) {
    mpi_parallel::cout << "+--------------------------------+" << std::endl;
    mpi_parallel::cout << "! Two level MC                   !" << std::endl;
    mpi_parallel::cout << "+--------------------------------+" << std::endl;
    mpi_parallel::cout << std::endl;

    std::shared_ptr<SamplerFactory> sampler_factory;
    sampler_factory = construct_sampler_factory(
        param_twolevelmc.sampler(), qoi_factory, coarse_sampler_factory,
        conditioned_fine_action_factory, param_general, param_qft, param_hmc,
        param_cluster, param_heatbath, param_hierarchical, param_stats);
    MonteCarloTwoLevel montecarlo_twolevel(action, qoi_factory, sampler_factory,
                                           conditioned_fine_action_factory,
                                           param_stats, param_twolevelmc);
    montecarlo_twolevel.evaluate_difference();
    montecarlo_twolevel.show_statistics();
    mpi_parallel::cout << std::endl;
  }

  /* **************************************** *
   * Multilevel method                         *
   * **************************************** */
  if (param_general.method() == MethodMultiLevel) {
    if (param_qft.action() == ActionNonlinearSigma) {
      mpi_parallel::cerr << " ERROR: multilevel method not yet supported for "
                            "nonlinear sigma model."
                         << std::endl;
      mpi_exit(EXIT_FAILURE);
    }
    if (mpi_comm_size() > 1) {
      mpi_parallel::cerr
          << " Multilevel method has not been parallelised (yet)." << std::endl;
      mpi_exit(EXIT_FAILURE);
    }
    mpi_parallel::cout << "+--------------------------------+" << std::endl;
    mpi_parallel::cout << "! Multilevel MC                  !" << std::endl;
    mpi_parallel::cout << "+--------------------------------+" << std::endl;
    mpi_parallel::cout << std::endl;

    std::shared_ptr<SamplerFactory> sampler_factory;
    sampler_factory = construct_sampler_factory(
        param_twolevelmc.sampler(), qoi_factory, coarse_sampler_factory,
        conditioned_fine_action_factory, param_general, param_qft, param_hmc,
        param_cluster, param_heatbath, param_hierarchical, param_stats);

    MonteCarloMultiLevel montecarlo_multilevel(
        action, qoi_factory, sampler_factory, conditioned_fine_action_factory,
        param_stats, param_multilevelmc);

    montecarlo_multilevel.evaluate();
    montecarlo_multilevel.show_statistics();
    if (param_multilevelmc.show_detailed_stats()) {
      montecarlo_multilevel.show_detailed_statistics();
    }
    numerical_result = montecarlo_multilevel.numerical_result();
    statistical_error = montecarlo_multilevel.statistical_error();
  }

  // print out comparison to analytical result, if this exists
  if (do_analytical_comparison) {
    if ((param_general.method() == MethodSingleLevel) or
        (param_general.method() == MethodMultiLevel)) {
      double diff = fabs(numerical_result - analytical_result);
      double ratio = diff / statistical_error;
      mpi_parallel::cout << std::setprecision(8) << std::fixed;
      mpi_parallel::cout << "Comparison to analytical result " << std::endl;
      mpi_parallel::cout << "  (analytical - numerical) = " << diff;
      mpi_parallel::cout << std::setprecision(3) << std::fixed;
      mpi_parallel::cout << " = " << ratio << " * (statistical error) "
                         << std::endl
                         << std::endl;
    }
  }
  total_time.stop();
  mpi_parallel::cout << total_time << std::endl;
  mpi_finalize();
}
