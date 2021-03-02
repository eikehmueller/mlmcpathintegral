#include <iostream>
#include <iomanip>
#include <utility>
#include <memory>

#include "common/parameters.hh"
#include "common/statistics.hh"
#include "mpi/mpi_wrapper.hh"
#include "lattice/lattice2d.hh"
#include "action/renormalisation.hh"
#include "action/qft/quenchedschwingeraction.hh"
#include "action/qft/quenchedschwingerconditionedfineaction.hh"
#include "qoi/qft/qoiavgplaquette.hh"
#include "qoi/qft/qoi2dsusceptibility.hh"
#include "montecarlo/montecarlosinglelevel.hh"
#include "montecarlo/montecarlotwolevel.hh"
#include "montecarlo/montecarlomultilevel.hh"
#include "sampler/hmcsampler.hh"
#include "sampler/hierarchicalsampler.hh"
#include "sampler/multilevelsampler.hh"
#include "sampler/overrelaxedheatbathsampler.hh"
#include "config.h"

/** @file driver_schwinger.cc
 * @brief File with main program for 2D Schwinger model
 *
 * @mainpage
 * Several classes for implementating Multilevel MCMC for the path-integral
 * formulation of the lattice Schwinger model in 2D.
 */

/** Helper function to construct suitable sampler factory for given samplerid */
std::shared_ptr<SamplerFactory> construct_sampler_factory(const int samplerid,
                                                          const std::shared_ptr<QoI2DSusceptibilityFactory> qoi_factory,
                                                          const std::shared_ptr<SamplerFactory> coarse_sampler_factory,
                                                          const std::shared_ptr<ConditionedFineActionFactory> conditioned_fine_action_factory,
                                                          const GeneralParameters param_general,
                                                          const HMCParameters param_hmc,
                                                          const OverrelaxedHeatBathParameters param_heatbath,
                                                          const HierarchicalParameters param_hierarchical,
                                                          const StatisticsParameters param_stats) {
    std::shared_ptr<SamplerFactory> sampler_factory;
    if (samplerid == SamplerHMC) {
        /* --- CASE 1: HMC sampler ---- */
        sampler_factory = std::make_shared<HMCSamplerFactory>(param_hmc);
    } else if (samplerid == SamplerOverrelaxedHeatBath) {
        /* --- CASE 2: heat bath sampler ---- */
        sampler_factory = std::make_shared<OverrelaxedHeatBathSamplerFactory>(param_heatbath);
    } else if (samplerid == SamplerHierarchical) {
        /* --- CASE 3: Hierarchical sampler */
        sampler_factory = std::make_shared<HierarchicalSamplerFactory>(coarse_sampler_factory,
                                                                       conditioned_fine_action_factory,
                                                                       param_hierarchical);
    } else if (samplerid == SamplerMultilevel) {
        /* --- CASE 4: Multilevel sampler */
        sampler_factory = std::make_shared<MultilevelSamplerFactory>(qoi_factory,
                                                                     coarse_sampler_factory,
                                                                     conditioned_fine_action_factory,
                                                                     param_stats,
                                                                     param_hierarchical);
    } else {
        mpi_parallel::cerr << " ERROR: Unsupported sampler." << std::endl;
        mpi_exit(EXIT_FAILURE);
    }
    return sampler_factory;
}


/** Main program */
int main(int argc, char* argv[]) {
    mpi_init();
    Timer total_time("total");
    total_time.start();
    mpi_parallel::cout << "++===================================++" << std::endl;
    mpi_parallel::cout << "!!   Path integral multilevel MCMC   !!" << std::endl;
    mpi_parallel::cout << "!!   for the 2D Schwinger model      !!" << std::endl;
    mpi_parallel::cout << "++===================================++" << std::endl;
    mpi_parallel::cout << std::endl;
#ifdef USE_MPI
    mpi_parallel::cout << "MPI parallel version running on " << mpi_comm_size() << " processes." << std::endl;
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
    mpi_parallel::cout << " Reading parameter from file \'" << filename << "\'" << std::endl;
    mpi_parallel::cout << std::endl;

    /* ====== Read parameters ====== */
    GeneralParameters param_general;
    Lattice2DParameters param_lattice;
    StatisticsParameters param_stats;
    SchwingerParameters param_schwinger;
    if (param_general.readFile(filename)) return 1;
    mpi_parallel::cout << param_general << std::endl;
    if (param_lattice.readFile(filename)) return 1;
    mpi_parallel::cout << param_lattice << std::endl;
    if (param_stats.readFile(filename)) return 1;
    mpi_parallel::cout << param_stats << std::endl;
    if (param_schwinger.readFile(filename)) return 1;
    mpi_parallel::cout << param_schwinger << std::endl;
    
    HMCParameters param_hmc;
    if (param_hmc.readFile(filename)) return 1;
    mpi_parallel::cout << param_hmc << std::endl;
  
    OverrelaxedHeatBathParameters param_heatbath;
    if (param_heatbath.readFile(filename)) return 1;
    mpi_parallel::cout << param_heatbath << std::endl;

    SingleLevelMCParameters param_singlelevelmc;
    if (param_singlelevelmc.readFile(filename)) return 1;
    mpi_parallel::cout << param_singlelevelmc << std::endl;
    
    HierarchicalParameters param_hierarchical;
    if (param_hierarchical.readFile(filename)) return 1;
    mpi_parallel::cout << param_hierarchical << std::endl;
    
    TwoLevelMCParameters param_twolevelmc;
    if (param_twolevelmc.readFile(filename)) return 1;
    mpi_parallel::cout << param_twolevelmc << std::endl;

    MultiLevelMCParameters param_multilevelmc;
    if (param_multilevelmc.readFile(filename)) return 1;
    mpi_parallel::cout << param_multilevelmc << std::endl;

#ifdef DEBUG_BUILD
    mpi_parallel::cout << FRED("CAUTION: built in debug mode.") << std::endl;
#endif // DEBUG_BUILD

#ifdef OPT_BUILD
    mpi_parallel::cout << FGREEN("Built in optimised mode.") << std::endl;
#endif // OPT_BUILD

#ifdef SAVE_PATHS
    mpi_parallel::cout << FRED("CAUTION: logging paths will impact performance!") << std::endl;
#endif // SAVE_PATHS

#ifdef LOG_QOI
    mpi_parallel::cout << FRED("CAUTION: logging QoI will impact performance!") << std::endl;
#endif // LOG_QOI

    /* ====== Lattice ====== */
    std::shared_ptr<Lattice2D> lattice;
    lattice = std::make_shared<Lattice2D>(param_lattice.Mt_lat(),
                                          param_lattice.Mx_lat(),
                                          param_lattice.T_lat(),
                                          param_lattice.L_lat());

    /* ====== Select quantity of interest ====== */
    std::shared_ptr<QoI2DSusceptibility> qoi;
    qoi = std::make_shared<QoI2DSusceptibility>(lattice);
    std::shared_ptr<QoI2DSusceptibilityFactory> qoi_factory = std::make_shared<QoI2DSusceptibilityFactory>();
    mpi_parallel::cout << std::endl;

    /* ====== Select action ====== */
    std::shared_ptr<Action> action;
    action = std::make_shared<QuenchedSchwingerAction>(lattice,
                                                       param_schwinger.renormalisation(),
                                                       param_schwinger.beta());

    // numerical result and statistical error
    double numerical_result;
    double statistical_error;
    double analytical_result = quenchedschwinger_chit_analytical(param_schwinger.beta(),
                                                                 lattice->getNcells());
    double analytical_result_variance = quenchedschwinger_var_chit_continuum_analytical(param_schwinger.beta(),
                                                                                        lattice->getNcells());
    mpi_parallel::cout << std::endl;
    mpi_parallel::cout << std::setprecision(8) << std::fixed;
    mpi_parallel::cout << " Analytical results"  << std::endl;
    mpi_parallel::cout << "      E[V*chi_t]              = " << analytical_result << std::endl;
    mpi_parallel::cout << "      lim_{a->0} Var[V*chi_t] = " << analytical_result_variance << std::endl;
    mpi_parallel::cout << std::endl;

    /* Construction conditioned fine action factory */
    std::shared_ptr<ConditionedFineActionFactory> conditioned_fine_action_factory
        = std::make_shared<QuenchedSchwingerConditionedFineActionFactory>();

    /* Construct coarse level sampler factory, which might be used by the hierarchical samplers */
    std::shared_ptr<SamplerFactory> coarse_sampler_factory;
    /* Note that here it does not make sense to use the hierarchical- or multilevel-sampler,
     * so we can pass null pointers for the QoI and coarse level sampler factory
     */
    coarse_sampler_factory = construct_sampler_factory(param_hierarchical.coarsesampler(),
                                                       nullptr,
                                                       nullptr,
                                                       nullptr,
                                                       param_general,
                                                       param_hmc,
                                                       param_heatbath,
                                                       param_hierarchical,
                                                       param_stats);
    
    /* **************************************** *
     * Single level method                      *
     * **************************************** */

    if (param_general.method() == MethodSingleLevel) {
        mpi_parallel::cout << "+--------------------------------+" << std::endl;
        mpi_parallel::cout << "! Single level MC                !" << std::endl;
        mpi_parallel::cout << "+--------------------------------+" << std::endl;
        mpi_parallel::cout << std::endl;

        std::shared_ptr<SamplerFactory> sampler_factory;
        sampler_factory = construct_sampler_factory(param_singlelevelmc.sampler(),
                                                    qoi_factory,
                                                    coarse_sampler_factory,
                                                    conditioned_fine_action_factory,
                                                    param_general,
                                                    param_hmc,
                                                    param_heatbath,
                                                    param_hierarchical,
                                                    param_stats);

        /* ====== Construct single level MC ====== */
        MonteCarloSingleLevel montecarlo_singlelevel(action,
                                                     qoi,
                                                     sampler_factory,
                                                     param_stats,
                                                     param_singlelevelmc);

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
        sampler_factory = construct_sampler_factory(param_twolevelmc.sampler(),
                                                    qoi_factory,
                                                    coarse_sampler_factory,
                                                    conditioned_fine_action_factory,
                                                    param_general,
                                                    param_hmc,
                                                    param_heatbath,
                                                    param_hierarchical,
                                                    param_stats);
        MonteCarloTwoLevel montecarlo_twolevel(action,
                                               qoi_factory,
                                               sampler_factory,
                                               conditioned_fine_action_factory,
                                               param_twolevelmc);
        Statistics stats_fine("QoI[fine]",10);
        Statistics stats_coarse("QoI[coarse]",10);
        Statistics stats_diff("delta QoI",10);
        montecarlo_twolevel.evaluate_difference(stats_fine,
                                                stats_coarse,
                                                stats_diff);
        mpi_parallel::cout << stats_fine << std::endl;
        mpi_parallel::cout << stats_coarse << std::endl;
        mpi_parallel::cout << stats_diff << std::endl;
        mpi_parallel::cout << std::endl;
        mpi_parallel::cout << "=== Coarse level sampler statistics === " << std::endl;
        montecarlo_twolevel.get_coarsesampler()->show_stats();
        mpi_parallel::cout << std::endl;
        mpi_parallel::cout << "=== Two level sampler statistics === " << std::endl;
        montecarlo_twolevel.get_twolevelstep()->show_stats();
        mpi_parallel::cout << std::endl;
    }
    
    /* **************************************** *
     * Multilevel method                         *
     * **************************************** */
    if (param_general.method() == MethodMultiLevel) {
        if (mpi_comm_size() > 1) {
            mpi_parallel::cerr << " Multilevel method has not been parallelised (yet)." << std::endl;
            mpi_exit(EXIT_FAILURE);
        }
        mpi_parallel::cout << "+--------------------------------+" << std::endl;
        mpi_parallel::cout << "! Multilevel MC                  !" << std::endl;
        mpi_parallel::cout << "+--------------------------------+" << std::endl;
        mpi_parallel::cout << std::endl;

        std::shared_ptr<SamplerFactory> sampler_factory;
        sampler_factory = construct_sampler_factory(param_twolevelmc.sampler(),
                                                    qoi_factory,
                                                    coarse_sampler_factory,
                                                    conditioned_fine_action_factory,
                                                    param_general,
                                                    param_hmc,
                                                    param_heatbath,
                                                    param_hierarchical,
                                                    param_stats);

        MonteCarloMultiLevel montecarlo_multilevel(action,
                                                   qoi_factory,
                                                   sampler_factory,
                                                   conditioned_fine_action_factory,
                                                   param_stats,
                                                   param_multilevelmc);

        montecarlo_multilevel.evaluate();
        montecarlo_multilevel.show_statistics();
        if (param_multilevelmc.show_detailed_stats()) {
            montecarlo_multilevel.show_detailed_statistics();
        }
        numerical_result = montecarlo_multilevel.numerical_result();
        statistical_error = montecarlo_multilevel.statistical_error();
    }
    
    if ( (param_general.method() == MethodSingleLevel) or
         (param_general.method() == MethodMultiLevel) ) {
        double diff = fabs(numerical_result-analytical_result);
        double ratio = diff/statistical_error;
        mpi_parallel::cout << std::setprecision(8) << std::fixed;
        mpi_parallel::cout << "Comparison to analytical result " << std::endl;
        mpi_parallel::cout << "  (analytical - numerical) = " << diff;
        mpi_parallel::cout << std::setprecision(3) << std::fixed;
        mpi_parallel::cout << " = " << ratio << " * (statistical error) " << std::endl << std::endl;
    }
    total_time.stop();
    mpi_parallel::cout << total_time << std::endl;
    mpi_finalize();
}
