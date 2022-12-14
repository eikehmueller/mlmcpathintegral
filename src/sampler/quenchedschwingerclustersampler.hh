#ifndef QUENCHEDSCHWINGERCLUSTERSAMPLER_HH
#define QUENCHEDSCHWINGERCLUSTERSAMPLER_HH QUENCHEDSCHWINGERCLUSTERSAMPLER_HH

#include "action/qft/quenchedschwingeraction.hh"
#include "action/qm/rotoraction.hh"
#include "action/renormalisation.hh"
#include "common/auxilliary.hh"
#include "common/samplestate.hh"
#include "common/timer.hh"
#include "lattice/lattice1d.hh"
#include "lattice/lattice2d.hh"
#include "mpi/mpi_random.hh"
#include "sampler/clustersampler.hh"
#include "sampler/sampler.hh"
#include <memory>
#include <random>

/** @file quenchedschwingerclustersampler.hh
 * @brief Header file for cluster sampler for quenched schwinger model
 */

/** @class QuenchedSchwingerClusterSampler
 * @brief Cluster sampler algorithm for the quenched Schwinger model
 *
 * The sampler exploits the equivalence between the quenched Schwinger model and
 * the topological oscillator to generate samples from the former with the
 * cluster sampler.
 *
 * The key observation is that the distribution of the differences
 * \f$x_j-x_{j-1}\f$ for a topological oscillator discretised with \f$M\f$ sites
 * is the same as the distribution of the plaquettes of the Schwinger model on a
 * lattice with \f$M = M_t\cdot M_x\f$ cells.
 *
 * This sampler uses the ClusterSampler class to generator samples from the
 * topological oscillator distribution, converts them to plaquettes which can be
 * mapped to configurations of links on the lattice.
 */
class QuenchedSchwingerClusterSampler : public Sampler {
public:
  /** @brief Create new instance
   *
   * @param[in] action_ Action to sample from
   * @param[in] cluster_param_ Parameters for the cluster sampler
   */
  QuenchedSchwingerClusterSampler(
      const std::shared_ptr<QuenchedSchwingerAction> action_,
      const ClusterParameters cluster_param);

  /** @brief Destroy instance
   *
   * Deallocate memory
   */
  virtual ~QuenchedSchwingerClusterSampler() {}

  /** @brief Draw a sample
   *
   * returns a sample state \f$\phi\f$
   *
   * @param[out] phi_state State \f$\phi\f$ drawn from distribution
   */
  virtual void draw(std::shared_ptr<SampleState> phi_state);

  /** Return cost per sample */
  virtual double cost_per_sample() { return cost_per_sample_; }

private:
  /** @brief Set current state to particular value
   *
   * @param[in] phi_state
   */
  virtual void set_state(std::shared_ptr<SampleState> phi_state){};

protected:
  /** @brief Action to sample from */
  const std::shared_ptr<QuenchedSchwingerAction> action;
  /** @brief State in 1d cluster sampler */
  mutable std::shared_ptr<SampleState> psi_cluster_state;
  /** @brief Underlying topological oscillator action */
  mutable std::shared_ptr<RotorAction> rotor_action;
  /** @brief underlying cluster sampler */
  mutable std::shared_ptr<ClusterSampler> cluster_sampler;
  /** @brief cost per sample */
  double cost_per_sample_;
  /** @brief Random number engine */
  mpi_parallel::mt19937_64 engine;
  /** @brief Uniform distribution for gauge transformation */
  std::uniform_real_distribution<double> uniform_distribution;
};

class QuenchedSchwingerClusterSamplerFactory : public SamplerFactory {
public:
  /** @brief Create new instance
   *
   * @param[in] param_cluster Custer sampler parameters
   */
  QuenchedSchwingerClusterSamplerFactory(const ClusterParameters param_cluster_)
      : param_cluster(param_cluster_) {}

  /** @brief Destructor */
  virtual ~QuenchedSchwingerClusterSamplerFactory() {}

  /** @brief Return sampler for a specific  action
   *
   * @param[in] action Action to sample from
   */
  virtual std::shared_ptr<Sampler> get(std::shared_ptr<Action> action) {
    return std::make_shared<QuenchedSchwingerClusterSampler>(
        std::dynamic_pointer_cast<QuenchedSchwingerAction>(action),
        param_cluster);
  }

private:
  /** Cluster sampler parameters */
  const ClusterParameters param_cluster;
};

#endif // QUENCHEDSCHWINGERCLUSTERSAMPLER_HH
