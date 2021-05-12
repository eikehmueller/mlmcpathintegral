#ifndef MONTECARLOTWOLEVEL_HH
#define MONTECARLOTWOLEVEL_HH MONTECARLOTWOLEVEL_HH
#include <utility>
#include <cmath>
#include <sstream>
#include <iostream>
#include <iomanip>
#include "config.h"
#include "common/statistics.hh"
#include "common/parameters.hh"
#include "mpi/mpi_wrapper.hh"
#include "sampler/sampler.hh"
#include "action/action.hh"
#include "action/conditionedfineaction.hh"
#include "qoi/quantityofinterest.hh"
#include "montecarlo/twolevelmetropolisstep.hh"
#include "montecarlo/montecarlo.hh"

/** @file montecarlotwolevel.hh
 * @brief Header file for two level Monte Carlo classes
 */

/** @class TwoLevelMCParameters
 *
 * @brief Class for storing parameters of two level Monte Carlo integrator.
 */
class TwoLevelMCParameters : public Parameters {
public:
    /** @brief Construct a new instance */
    TwoLevelMCParameters() :
        Parameters("twolevelmc"),
        n_burnin_(100),
        n_samples_(100),
        n_coarse_autocorr_window_(10),
        n_fine_autocorr_window_(10),
        n_delta_autocorr_window_(10),
        sampler_(SamplerHMC) {
        addKey("n_burnin",Integer,Positive);
        addKey("n_samples",Integer,Positive);
        addKey("n_coarse_autocorr_window",Integer,Positive);
        addKey("n_fine_autocorr_window",Integer,Positive);
        addKey("n_delta_autocorr_window",Integer,Positive);
        addKey("sampler",String);
    }

    /** @brief Read parameters from file
     *
     * @param[in] filename Name of file to read
     */
    int readFile(const std::string filename) {

        int readSuccess = Parameters::readFile(filename);
        if (!readSuccess) {
            n_burnin_ = getContents("n_burnin")->getInt();
            n_samples_ = getContents("n_samples")->getInt();
            n_coarse_autocorr_window_ = getContents("n_coarse_autocorr_window")->getInt();
            n_fine_autocorr_window_ = getContents("n_fine_autocorr_window")->getInt();
            n_delta_autocorr_window_ = getContents("n_delta_autocorr_window")->getInt();
            std::string sampler_str = getContents("sampler")->getString();
            if (sampler_str == "HMC") {
                sampler_ = SamplerHMC;
            } else if (sampler_str == "heatbath") {
                sampler_ = SamplerOverrelaxedHeatBath;
            } else if (sampler_str == "cluster") {
                sampler_ = SamplerCluster;
            } else if (sampler_str == "hierarchical") {
                sampler_ = SamplerHierarchical;
            } else if (sampler_str == "exact") {
                sampler_ = SamplerExact;
            } else  {
                mpi_parallel::cerr << " ERROR: Unknown sampler: " << sampler_str;
                mpi_parallel::cerr << std::endl;
                mpi_parallel::cerr << "        allowed values are \'HMC\', \'heatbath\', \'cluster\', \'hierarchical\', \'exact\'" << std::endl;
                mpi_exit(EXIT_FAILURE);
            }
        }
        return readSuccess;
    }

    /** @brief Return number of burnin samples */
    unsigned int n_burnin() const {
        return n_burnin_;
    }
    /** @brief Return number of samples param_twolevelmc.n_autocorr_window()*/
    unsigned int n_samples() const {
        return n_samples_;
    }
    /** @brief Return size of autocorrelatin window for coarse level sampler */
    unsigned int n_coarse_autocorr_window() const {
        return n_coarse_autocorr_window_;
    }
    /** @brief Return size of autocorrelatin window for fine level */
    unsigned int n_fine_autocorr_window() const {
        return n_fine_autocorr_window_;
    }
    /** @brief Return size of autocorrelatin window for difference */
    unsigned int n_delta_autocorr_window() const {
        return n_delta_autocorr_window_;
    }
    /** @brief Return sampler type */
    SamplerType sampler() const {
        return sampler_;
    }
private:
    /** @brief Number of burnin samples */
    unsigned int n_burnin_;
    /** @brief Number of samples */
    unsigned int n_samples_;
    /** @brief Size of window used for measuring coarse level sampler autocorrelations */
    unsigned int n_coarse_autocorr_window_;
    /** @brief Size of window used for measuring fine level autocorrelations */
    unsigned int n_fine_autocorr_window_;
    /** @brief Size of window used for measuring difference autocorrelations */
    unsigned int n_delta_autocorr_window_;
    /** @brief Sampler type */
    SamplerType sampler_;
};

/** @class MonteCarloTwoLevel
 *
 * @brief Two level Monte Carlo method
 *
 *
 */
class MonteCarloTwoLevel : public MonteCarlo {
public:
    /** \brief Create new instance
     *
     * @param[in] fine_action_ Action on fine level
     * @param[in] qoi_ Quantity of interest
     * @param[in] sampler_factory Factory for constructing coarse level sampler
     * @param[in] sampler_factory Factory for constructing conditioned fine actions
     * @param[in] param_twolevelmc Two level sampler parameters
     */
    MonteCarloTwoLevel(const std::shared_ptr<Action> fine_action_,
                       const std::shared_ptr<QoIFactory> qoi_factory_,
                       const std::shared_ptr<SamplerFactory> sampler_factory,
                       const std::shared_ptr<ConditionedFineActionFactory> conditioned_fine_action_factory,
                       const StatisticsParameters param_stats,
                       const TwoLevelMCParameters param_twolevelmc);

    /** @brief Calculate mean and variance of difference
     *
     * Calculate the mean and variance of the difference in the QoI
     * evaluated at two subsequent levels.
     */
    void evaluate_difference();

    /** @brief Show statistics */
    void show_statistics() const;

private:
    
    /** @brief Draw independent sample on coarse level
     *
     * By subsampling the coarse level sampler, draw a new approximately
     * independent coarse level sample \f$\phi^{(c)}\f$
     *
     * @param[inout] phi_state Coarse level state \f$\phi^{(c)}\f$
     */
    void draw_coarse_sample(std::shared_ptr<SampleState> phi_state);
    
    /** @brief Number of samples */
    const unsigned int n_samples;
    /** @brief Sampler on coarse level */
    std::shared_ptr<Sampler> coarse_sampler;
    /** @brief Action on coarse level */
    std::shared_ptr<Action> coarse_action;
    /** @brief Action on fine level */
    std::shared_ptr<Action> fine_action;
    /** @brief Conditioned fine action */
    std::shared_ptr<ConditionedFineAction> conditioned_fine_action;
    /** @brief Quantity of interest on fine level */
    std::shared_ptr<QoI> qoi_fine;
    /** @brief Quantity of interest on coarse level */
    std::shared_ptr<QoI> qoi_coarse;
    /** Two-level sampler */
    std::shared_ptr<TwoLevelMetropolisStep> twolevel_step;
    /** @brief number of skipped samples between independent coarse samples on coarse levels */
    double t_indep;
    /** @brief number of independent coarse samples on coarse level */
    int n_indep;
    /** @brief number of samples generated since last independent coarse sample */
    int t_sampler;
    /** @brief Statistics on fine level */
    Statistics stats_fine;
    /** @brief Statistics on coarse level */
    Statistics stats_coarse;
    /** @brief Statistics on coarse level */
    Statistics stats_diff;
    /** @brief Statistics of coarse sampler */
    Statistics stats_coarse_sampler;
};

#endif // MONTECARLOTWOLEVEL_HH
