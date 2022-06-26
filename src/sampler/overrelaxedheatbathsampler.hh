#ifndef OVERRELAXEDHEATBATHSAMPLER_HH
#define OVERRELAXEDHEATBATHSAMPLER_HH OVERRELAXEDHEATBATHSAMPLER_HH
#include "action/action.hh"
#include "common/parameters.hh"
#include "common/samplestate.hh"
#include "mpi/mpi_random.hh"
#include "mpi/mpi_wrapper.hh"
#include "sampler/sampler.hh"
#include <memory>
#include <numeric>
#include <random>
#include <vector>

/** @file overrelaxedheatbathsampler.hh
 * @brief Header file for overrelaxed heat bath sampler class
 */

/** @class HeatBathParameters
 *
 * @brief Class for storing parameters of heat bath sampler
 */
class OverrelaxedHeatBathParameters : public Parameters {
public:
  /** @brief Construct a new instance */
  OverrelaxedHeatBathParameters()
      : Parameters("heatbath"), n_sweep_heatbath_(1), n_sweep_overrelax_(1),
        n_burnin_(100), random_order_(true) {
    addKey("n_sweep_heatbath", Integer, Positive);
    addKey("n_sweep_overrelax", Integer, Positive);
    addKey("n_burnin", Integer, Positive);
    addKey("random_order", Bool);
  }

  /** @brief Read parameters from file
   *
   * @param[in] filename Name of file to read
   */
  int readFile(const std::string filename) {

    int readSuccess = Parameters::readFile(filename);
    if (!readSuccess) {
      n_sweep_heatbath_ = getContents("n_sweep_heatbath")->getInt();
      n_sweep_overrelax_ = getContents("n_sweep_overrelax")->getInt();
      n_burnin_ = getContents("n_burnin")->getInt();
      random_order_ = getContents("random_order")->getBool();
    }
    return readSuccess;
  }

  /** @brief Return number of heatbath sweeps */
  unsigned int n_sweep_heatbath() const { return n_sweep_heatbath_; }

  /** @brief Return number of overrelaxation sweeps */
  unsigned int n_sweep_overrelax() const { return n_sweep_overrelax_; }

  /** @brief Return number of burnin samples */
  unsigned int n_burnin() const { return n_burnin_; }

  /** @brief Return true if unknowns are traversed in random order */
  bool random_order() const { return random_order_; }

private:
  /** @brief Number of heat bath sweeps */
  unsigned int n_sweep_heatbath_;
  /** @brief Number of overrelaxation sweeps */
  unsigned int n_sweep_overrelax_;
  /** @brief Number of burnin samples */
  unsigned int n_burnin_;
  /** @brief Traverse unknowns randomly in a sweep? */
  bool random_order_;
};

/** @class OverrelaxedHeatBathSampler
 * @brief Overrelaxed heat bath sampler for gauge theories.
 *
 * Implementation of heat bath sampling with overrelaxation for gauge theories.
 * Generate states by drawing from the local action, conditioned on the
 * neighbouring sites. For this, pass over the unknowns n_{sweep,overrelax}
 * times and apply the overrelaxation method at each site. Following this, carry
 * out n_{sweep,heatbath} sweeps over the unknownes, and locally upate each site
 * with a heat bath step.
 *
 * Sweeps over the unknowns can be carried out lexicographically or in random
 * order.
 *
 * References:
 *  * Creutz, M., 1987. "Monte Carlo algorithms for lattice gauge theory" (No.
 * BNL-39747). Brookhaven National Lab.
 *  * Brown, F.R. and Woch, T.J., 1987. "Overrelaxed heat-bath and Metropolis
 * algorithms for accelerating pure gauge Monte Carlo calculations." Physical
 * Review Letters, 58(23), p.2394.
 *  * Degrand, T.A. and DeTar, C., 2006. "Lattice methods for quantum
 * chromodynamics." World Scientific (Chapter 7)
 */
class OverrelaxedHeatBathSampler : public Sampler {
public:
  /** @brief Create new instance
   *
   * @param[in] action_ Action to sample from
   * @param[in] heatbath_param_ Heat bath parameters
   */
  OverrelaxedHeatBathSampler(
      const std::shared_ptr<Action> action_,
      const OverrelaxedHeatBathParameters heatbath_param_)
      : Sampler(), action(action_),
        n_sweep_heatbath(heatbath_param_.n_sweep_heatbath()),
        n_sweep_overrelax(heatbath_param_.n_sweep_overrelax()),
        n_burnin(heatbath_param_.n_burnin()),
        random_order(heatbath_param_.random_order()),
        index_map(action->get_heatbath_indexset()) {
    engine.seed(871417);
    // If the index map is empty, fill it with consecutive numbers
    if (index_map.empty()) {
      index_map.resize(action->sample_size());
      std::iota(index_map.begin(), index_map.end(), 0);
    }
    // Create temporary workspace
    // current state
    phi_state_cur = std::make_shared<SampleState>(action->sample_size());
    action->initialise_state(phi_state_cur);
    // Burn in
    std::shared_ptr<SampleState> phi_state_tmp =
        std::make_shared<SampleState>(action->sample_size());
    for (unsigned int i = 0; i < n_burnin; ++i) {
      draw(phi_state_tmp);
    }
    // Reset counters
    reset_stats();
  }

  /** @brief Destroy instance
   *
   * Deallocate memory
   */
  virtual ~OverrelaxedHeatBathSampler() {}

  /** @brief Draw a new sample
   *
   * returns a sample state \f$\phi\f$
   *
   * @param[out] phi_state Path \f$\phi\f$ drawn from distribution
   */
  virtual void draw(std::shared_ptr<SampleState> phi_state);

  /** @brief Set current state to particular value
   *
   * @param[in] phi_state
   */
  virtual void set_state(std::shared_ptr<SampleState> phi_state);

protected:
  /** @brief Action to sample from */
  const std::shared_ptr<Action> action;
  /** @brief Number of heat bath sweeps over unknowns */
  const unsigned int n_sweep_heatbath;
  /** @brief Number of overrelaxation sweeps over unknowns */
  const unsigned int n_sweep_overrelax;
  /** @brief Number of burn-in steps */
  const unsigned int n_burnin;
  /** @brief Traverse unknowns randomly in a sweep? */
  bool random_order;
  /** @brief Current state */
  mutable std::shared_ptr<SampleState> phi_state_cur;
  /** @brief Type of Mersenne twister engine */
  typedef mpi_parallel::mt19937_64 Engine;
  /** @brief Random number engine used for shuffling (if unknowns are traversed
   * in random order) */
  mutable Engine engine;
  /** @brief Indirection map, used for lattice traversal in random order */
  std::vector<unsigned int> index_map;
};

class OverrelaxedHeatBathSamplerFactory : public SamplerFactory {
public:
  /** @brief Create new instance
   *
   * @param[in] param_heatbath Heat bath parameters
   */
  OverrelaxedHeatBathSamplerFactory(
      const OverrelaxedHeatBathParameters param_heatbath_)
      : param_heatbath(param_heatbath_) {}

  /** @brief Destructor */
  virtual ~OverrelaxedHeatBathSamplerFactory() {}

  /** @brief Return sampler for a specific  action
   *
   * @param[in] action Action to sample from
   */
  virtual std::shared_ptr<Sampler> get(std::shared_ptr<Action> action) {
    return std::make_shared<OverrelaxedHeatBathSampler>(action, param_heatbath);
  }

private:
  /** Heat bath parameters */
  const OverrelaxedHeatBathParameters param_heatbath;
};

#endif // OVERRELAXEDHEATBATHSAMPLER_HH
