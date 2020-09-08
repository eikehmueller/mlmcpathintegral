#include <iostream>
#include <iomanip>
#include <utility>
#include <memory>

#include "common/parameters.hh"
#include "common/statistics.hh"
#include "mpi/mpi_wrapper.hh"
#include "lattice/lattice1d.hh"
#include "action/renormalisation.hh"
#include "action/qm/qmaction.hh"
#include "action/qm/harmonicoscillatoraction.hh"
#include "action/qm/quarticoscillatoraction.hh"
#include "action/qm/rotoraction.hh"
#include "action/qm/gaussianconditionedfineaction.hh"
#include "action/qm/rotorconditionedfineaction.hh"
#include "qoi/qm/qoixsquared.hh"
#include "qoi/qm/qoisusceptibility.hh"
#include "montecarlo/montecarlosinglelevel.hh"
#include "montecarlo/montecarlotwolevel.hh"
#include "montecarlo/montecarlomultilevel.hh"
#include "sampler/hmcsampler.hh"
#include "sampler/clustersampler.hh"
#include "sampler/multilevelsampler.hh"
#include "config.h"

/** @file driver.cc
 * @brief File with main program
 *
 * @mainpage
 * Several classes for implementating Multilevel MCMC for the path-integral
 * formulation of quantum mechanics.
 */

/** Helper function to construct suitable sampler factory for given samplerid */
std::shared_ptr<SamplerFactory> construct_sampler_factory(const int samplerid,
        const std::shared_ptr<QoI> qoi,
        const std::shared_ptr<SamplerFactory> coarse_sampler_factory,
        const std::shared_ptr<ConditionedFineActionFactory> conditioned_fine_action_factory,
        const GeneralParameters param_general,
        const QMParameters param_qm,
        const HMCParameters param_hmc,
        const ClusterParameters param_cluster,
        const StatisticsParameters param_stats,
        const HierarchicalParameters param_hierarchical) {
    std::shared_ptr<SamplerFactory> sampler_factory;
    if (samplerid == SamplerHMC) {
        /* --- CASE 1: HMC sampler ---- */
        sampler_factory = std::make_shared<HMCSamplerFactory>(param_hmc);
    } else if (samplerid == SamplerCluster) {
        /* --- CASE 2: cluster sampler ---- */
        if (param_qm.action() != ActionRotor) {
            mpi_parallel::cerr << " ERROR: can only use cluster sampler for QM rotor action." << std::endl;
            mpi_exit(EXIT_FAILURE);
        }
        sampler_factory = std::make_shared<ClusterSamplerFactory>(param_cluster);
    } else if (samplerid == SamplerExact) {
        /* --- CASE 3: exact sampler (for HO action) ---- */
        if (param_qm.action() != ActionHarmonicOscillator) {
            mpi_parallel::cerr << " ERROR: can only sample exactly from harmonic oscillator action." << std::endl;
            mpi_exit(EXIT_FAILURE);
        }
        sampler_factory = std::make_shared<HarmonicOscillatorSamplerFactory>();
    } else if (samplerid == SamplerHierarchical) {
        /* ---- CASE 4: Hierarchical sampler */
        sampler_factory = std::make_shared<HierarchicalSamplerFactory>(coarse_sampler_factory,
                          conditioned_fine_action_factory,
                          param_hierarchical);
    } else if (samplerid == SamplerMultilevel) {
        /* ---- CASE 5: Multilevel sampler */
        sampler_factory = std::make_shared<MultilevelSamplerFactory>(qoi,
                          coarse_sampler_factory,
                          conditioned_fine_action_factory,
                          param_stats,
                          param_hierarchical);
    } else {
        mpi_parallel::cerr << " ERROR: Unknown sampler." << std::endl;
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
    QMParameters param_qm;
    LatticeParameters param_lattice;
    StatisticsParameters param_stats;
    HarmonicOscillatorParameters param_ho;
    QuarticOscillatorParameters param_qo;
    RotorParameters param_rotor;
    if (param_general.readFile(filename)) return 1;
    mpi_parallel::cout << param_general << std::endl;
    if (param_qm.readFile(filename)) return 1;
    mpi_parallel::cout << param_qm << std::endl;
    if (param_lattice.readFile(filename)) return 1;
    mpi_parallel::cout << param_lattice << std::endl;
    if (param_stats.readFile(filename)) return 1;
    mpi_parallel::cout << param_stats << std::endl;
    switch (param_qm.action()) {
    case (ActionHarmonicOscillator): {
        if (param_ho.readFile(filename)) return 1;
        mpi_parallel::cout << param_ho << std::endl;
        break;
    }
    case (ActionQuarticOscillator): {
        if (param_qo.readFile(filename)) return 1;
        mpi_parallel::cout << param_qo << std::endl;
        break;
    }
    case (ActionRotor): {
        if (param_rotor.readFile(filename)) return 1;
        mpi_parallel::cout << param_rotor << std::endl;
        break;
    }
    }

    HMCParameters param_hmc;
    if (param_hmc.readFile(filename)) return 1;
    mpi_parallel::cout << param_hmc << std::endl;

    ClusterParameters param_cluster;
    if (param_cluster.readFile(filename)) return 1;
    mpi_parallel::cout << param_cluster << std::endl;

    SingleLevelMCParameters param_singlelevelmc;
    if (param_singlelevelmc.readFile(filename)) return 1;
    mpi_parallel::cout << param_singlelevelmc << std::endl;

    TwoLevelMCParameters param_twolevelmc;
    if (param_twolevelmc.readFile(filename)) return 1;
    mpi_parallel::cout << param_twolevelmc << std::endl;

    MultiLevelMCParameters param_multilevelmc;
    if (param_multilevelmc.readFile(filename)) return 1;
    mpi_parallel::cout << param_multilevelmc << std::endl;

    HierarchicalParameters param_hierarchical;
    if (param_hierarchical.readFile(filename)) return 1;
    mpi_parallel::cout << param_hierarchical << std::endl;


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
    std::shared_ptr<Lattice1D> lattice;
    lattice = std::make_shared<Lattice1D>(param_lattice.M_lat(),
                                          param_lattice.T_final());

    /* ====== Select quantity of interest ====== */
    std::shared_ptr<QoI> qoi;
    mpi_parallel::cout << std::endl;
    if ( (param_qm.action() == ActionHarmonicOscillator) or
            (param_qm.action() == ActionQuarticOscillator) ) {
        qoi=std::make_shared<QoIXsquared>(lattice);
        mpi_parallel::cout << "QoI = X^2 " << std::endl;
    }
    if ( (param_qm.action() == ActionRotor) ) {
        qoi=std::make_shared<QoISusceptibility>(lattice);
        mpi_parallel::cout << "QoI = Susceptibility Q[X]^2/T " << std::endl;
    }
    mpi_parallel::cout << std::endl;

    /* ====== Select action ====== */
    std::shared_ptr<QMAction> action;
    switch (param_qm.action()) {
    case (ActionHarmonicOscillator): {
        action =
            std::make_shared<HarmonicOscillatorAction>(lattice,
                    param_ho.renormalisation(),
                    param_ho.m0(),
                    param_ho.mu2());
        break;
    }
    case (ActionQuarticOscillator): {
        action =
            std::make_shared<QuarticOscillatorAction>(lattice,
                    RenormalisationNone,
                    param_qo.m0(),
                    param_qo.mu2(),
                    param_qo.lambda(),
                    param_qo.x0());
        break;
    }
    case (ActionRotor): {
        action =
            std::make_shared<RotorAction>(lattice,
                                          param_rotor.renormalisation(),
                                          param_rotor.m0());

        break;
    }
    }

    // Exact result (if known)
    double exact_result;
    // numerical result and statistical error
    double numerical_result;
    double statistical_error;

    if ( (param_general.method() == MethodSingleLevel) or
            (param_general.method() == MethodMultiLevel) ) {
        /* ====== Print out exact result for harmonic oscillator */
        if (param_qm.action() == ActionHarmonicOscillator) {
            std::shared_ptr<HarmonicOscillatorAction> ho_action =
                std::dynamic_pointer_cast<HarmonicOscillatorAction>(action);
            exact_result = ho_action->Xsquared_exact();
            double exact_result_continuum = ho_action->Xsquared_exact_continuum();
            mpi_parallel::cout << std::endl;
            mpi_parallel::cout << std::setprecision(6) << std::fixed;
            mpi_parallel::cout << " Exact result             <x^2> = " << exact_result << std::endl;
            mpi_parallel::cout << " Continuum limit [a -> 0] <x^2> = " << exact_result_continuum << std::endl;
            mpi_parallel::cout << std::endl;
        }
        if (param_qm.action() == ActionRotor) {
            std::shared_ptr<RotorAction> rotor_action =
                std::dynamic_pointer_cast<RotorAction>(action);
            double exact_result_continuum = rotor_action->chit_exact_continuum();
            exact_result = rotor_action->chit_exact_perturbative();
            mpi_parallel::cout << std::endl;
            mpi_parallel::cout << std::setprecision(6) << std::fixed;
            mpi_parallel::cout << " Exact result             <chi_t> = " << exact_result << " + O((a/I)^2), a/I = " << lattice->geta_lat()/action->getm0() << std::endl;
            mpi_parallel::cout << " Continuum limit [a -> 0] <chi_t> = " << exact_result_continuum << std::endl;
            mpi_parallel::cout << std::endl;
        }
    }

    /* Construction conditioned fine action factory */
    std::shared_ptr<ConditionedFineActionFactory> conditioned_fine_action_factory;
    if (param_qm.action() == ActionRotor) {
        conditioned_fine_action_factory = std::make_shared<RotorConditionedFineActionFactory>();
    } else {
        conditioned_fine_action_factory = std::make_shared<GaussianConditionedFineActionFactory>();
    }

    /* Construct coarse level sampler factory, which might be used by the hierarchical samplers */
    std::shared_ptr<SamplerFactory> coarse_sampler_factory;
    /* Note that here it does not make sense to use the hierarchical- or multilevel-sampler,
     * so we can pass null pointers for the QoI and coarse level sampler factory
     */
    coarse_sampler_factory = construct_sampler_factory(param_hierarchical.coarsesampler(),
                             nullptr,
                             nullptr,
                             conditioned_fine_action_factory,
                             param_general,
                             param_qm,
                             param_hmc,
                             param_cluster,
                             param_stats,
                             param_hierarchical);

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
                          qoi,
                          coarse_sampler_factory,
                          conditioned_fine_action_factory,
                          param_general,
                          param_qm,
                          param_hmc,
                          param_cluster,
                          param_stats,
                          param_hierarchical);
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
                          qoi,
                          coarse_sampler_factory,
                          conditioned_fine_action_factory,
                          param_general,
                          param_qm,
                          param_hmc,
                          param_cluster,
                          param_stats,
                          param_hierarchical);
        MonteCarloTwoLevel montecarlo_twolevel(action,
                                               qoi,
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
        sampler_factory = construct_sampler_factory(param_multilevelmc.sampler(),
                          qoi,
                          coarse_sampler_factory,
                          conditioned_fine_action_factory,
                          param_general,
                          param_qm,
                          param_hmc,
                          param_cluster,
                          param_stats,
                          param_hierarchical);

        MonteCarloMultiLevel montecarlo_multilevel(action,
                qoi,
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

    // Compare numerical results to exact result (if possible)
    if ( ( (param_qm.action() == ActionHarmonicOscillator) or
            (param_qm.action() == ActionRotor) ) and
            ( (param_general.method() == MethodSingleLevel) or
              ( (param_general.method() == MethodMultiLevel) ) ) ) {
        double diff = fabs(numerical_result-exact_result);
        double ratio = diff/statistical_error;
        mpi_parallel::cout << std::setprecision(8) << std::fixed;
        mpi_parallel::cout << "Comparison to exact result " << std::endl;
        mpi_parallel::cout << "  (exact - numerical) = " << diff;
        mpi_parallel::cout << std::setprecision(3) << std::fixed;
        mpi_parallel::cout << " = " << ratio << " * (statistical error) " << std::endl << std::endl;
    }
    total_time.stop();
    mpi_parallel::cout << total_time << std::endl;
    mpi_finalize();
}
