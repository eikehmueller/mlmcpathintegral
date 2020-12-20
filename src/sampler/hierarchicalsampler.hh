#ifndef HIERARCHICALSAMPLER_HH
#define HIERARCHICALSAMPLER_HH HIERARCHICALSAMPLER_HH
#include <vector>
#include <memory>
#include <string>
#include <sstream>
#include "common/parameters.hh"
#include "common/statistics.hh"
#include "common/samplestate.hh"
#include "common/timer.hh"
#include "mpi/mpi_wrapper.hh"
#include "action/action.hh"
#include "sampler/sampler.hh"
#include "action/conditionedfineaction.hh"
#include "montecarlo/twolevelmetropolisstep.hh"

/** @class HierarchicalParameters
 *
 * @brief Class for storing parameters of Hierarchical sampler.
 */
class HierarchicalParameters : public Parameters {
public:
    /** @brief Construct a new instance */
    HierarchicalParameters() :
        Parameters("hierarchical"),
        n_maphi_level_(2),
        coarsesampler_(SamplerHMC) {
        addKey("n_maphi_level",Integer,Positive);
        addKey("coarsesampler",String);
    }

    /** @brief Read parameters from file
     *
     * @param[in] filename Name of file to read
     */
    int readFile(const std::string filename) {

        int readSuccess = Parameters::readFile(filename);
        if (!readSuccess) {
            n_maphi_level_ = getContents("n_maphi_level")->getInt();
            std::string sampler_str = getContents("coarsesampler")->getString();
            if (sampler_str == "HMC") {
                coarsesampler_ = SamplerHMC;
            } else if (sampler_str == "cluster") {
                coarsesampler_ = SamplerCluster;
            } else if (sampler_str == "exact") {
                coarsesampler_ = SamplerExact;
            } else  {
                mpi_parallel::cerr << " ERROR: Unknown coarse sampler: " << sampler_str;
                mpi_parallel::cerr << std::endl;
                mpi_parallel::cerr << "        allowed values are \'HMC\', \'cluster\', \'exact\'" << std::endl;
                mpi_exit(EXIT_FAILURE);
            }
        }
        return readSuccess;
    }

    /** @brief Return maximal number of levels */
    unsigned int n_max_level() const {
        return n_maphi_level_;
    }
    /** @brief Return sampler type */
    SamplerType coarsesampler() const {
        return coarsesampler_;
    }
private:
    /** @brief Number of levels */
    unsigned int n_maphi_level_;
    /** @brief tolerance epsilon */
    SamplerType coarsesampler_;
};

/** @file hierarchicalsampler.hh
 * @brief Header file for Hierarchical sampler class
 */

/** @class HierarchicalSampler
 * @brief Hierarchical sampler
 *
 * Implementation of hierarchical sampling. Generate samples by successively screening samplers
 * from coarser levels
 */
class HierarchicalSampler : public Sampler {
public:
    /** @brief Create new instance
       *
     * @param[in] fine_action_ Fine level Action to sample from
     * @param[in] coarse_sampler_factory Factory for generating coarse level sampler
     * @param[in] conditioned_fine_action_factory factory for generating conditioned fine actions
     * @param[in] param_hierarchical Hierarchical sampler parameters
     */
    HierarchicalSampler(const std::shared_ptr<Action> fine_action,
                        const std::shared_ptr<SamplerFactory> coarse_sampler_factory,
                        const std::shared_ptr<ConditionedFineActionFactory> conditioned_fine_action_factory,
                        const HierarchicalParameters param_hierarchical);

    /** @brief Destroy instance
     *
     * Deallocate memory
     */
    virtual ~HierarchicalSampler() {}

    /** @brief Draw a sample
     *
     * returns a sample state \f$\phi\f$
     *
     * @param[out] phi_state State \f$\phi\f$ drawn from distribution
     */
    virtual void draw(std::shared_ptr<SampleState> phi_state);

    /** @brief Show statistics
         *
         * Print out statistics
         */
    virtual void show_stats();

    /** @brief Set current state to particular value
     *
     * @param[in] phi_state
     */
    virtual void set_state(std::shared_ptr<SampleState> phi_state);

    /** Return cost per sample */
    virtual double cost_per_sample() {
        return cost_per_sample_;
    }

private:
    /** @brief Action to sample from */
    std::vector<std::shared_ptr<Action>> action;
    /** @brief Number of levels in hierarchy */
    const unsigned int n_level;
    /** @brief Two level step on all levels of the multigrid hierarchy */
    std::vector<std::shared_ptr<TwoLevelMetropolisStep> > twolevel_step;
    /** @brief Sampler on coarsest level */
    std::shared_ptr<Sampler> coarse_sampler;
    /** @brief Sampler state on each level  */
    std::vector<std::shared_ptr<SampleState> > phi_sampler_state;
    /** @brief cost per sample */
    double cost_per_sample_;
};

class HierarchicalSamplerFactory : public SamplerFactory {
public:
    /** @brief Create new instance
     *
     * @param[in] coarse_sampler_factory_ Factory for coarse level sampler
     * @param[in] conditioned_fine_action_factory factory for generating conditioned fine actions
     * @param[in] param_hierarchical Hierarchical sampler parameters
     */
    HierarchicalSamplerFactory(const std::shared_ptr<SamplerFactory> coarse_sampler_factory_,
                               const std::shared_ptr<ConditionedFineActionFactory> conditioned_fine_action_factory_,
                               const HierarchicalParameters param_hierarchical_) :
        coarse_sampler_factory(coarse_sampler_factory_),
        conditioned_fine_action_factory(conditioned_fine_action_factory_),
        param_hierarchical(param_hierarchical_) {}

    /** @brief Destructor */
    virtual ~HierarchicalSamplerFactory() {}

    /** @brief Return sampler for a specific  action
     *
     * @param[in] action Action to sample from
     */
    virtual std::shared_ptr<Sampler> get(std::shared_ptr<Action> action) {
        return std::make_shared<HierarchicalSampler>(action,
                coarse_sampler_factory,
                conditioned_fine_action_factory,
                param_hierarchical);
    }
private:
    /** Factory for coarsest level sampler*/
    const std::shared_ptr<SamplerFactory> coarse_sampler_factory;
    /** Conditioned fine action factory */
    const std::shared_ptr<ConditionedFineActionFactory> conditioned_fine_action_factory;
    /** Hierarchical sampler parameters */
    const HierarchicalParameters param_hierarchical;
};

#endif // HIERARCHICALSAMPLER
