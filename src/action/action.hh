#ifndef ACTION_HH
#define ACTION_HH ACTION_HH
#include "action/renormalisation.hh"
#include "common/samplestate.hh"
#include "mpi/mpi_wrapper.hh"
#include <cassert>
#include <iostream>
#include <memory>
#include <random>
#include <sstream>
#include <string>
#include <vector>

/** @file action.hh
 * @brief Header file for action base class
 */

/** @class Action
 *
 * @brief Base class for action
 *
 * Allows calculation of action
 *   \f[
        S[X]=\sum_{j=0}^{M-1}\left(\frac{m_0}{2}\frac{(X_{j+1}-X_j)}{a^2}+V(x)\right)
      \f]
 * for one-dimensional quantum problem with periodic boundary conditions.
 */
class Action {
public:
  /** @brief Initialise class
   *
   * @param[in] renormalisation_ Type of renormalisation
   */
  Action(const RenormalisationType renormalisation_)
      : renormalisation(renormalisation_) {}

  /** @brief return size of samples */
  virtual unsigned int sample_size() const = 0;

  /** @brief Cost of one action evaluation */
  virtual double evaluation_cost() const = 0;

  /** @brief Construct coarsened version of action
   *
   * This returns a coarsened version of the action on the next level
   * of the multigrid hierarchy.
   */
  virtual std::shared_ptr<Action> coarse_action() {
    mpi_parallel::cerr << "ERROR: cannot coarsen action" << std::endl;
    mpi_exit(EXIT_FAILURE);
    throw std::runtime_error("...");
  };

  /** @brief Evaluate action for a specific state
   *
   * Calculate \f$S[\phi]\f$ for a specific state
   *
   * @param[in] phi_state Sample state \f$\phi\f$
   */
  virtual const double
  evaluate(const std::shared_ptr<SampleState> phi_state) const = 0;

  /** @brief Draw local value of state from heat bath
   *
   * Update the local entry at position ell of the state using a local heatbath
   * (or Gibbs) step. Usually the index ell corresponds to a dof-index,
   * but depending on the action it can also be the index of a topological
   * entity.
   *
   *  @param[inout] phi_state State to update
   *  @param[in] ell index to update
   */
  virtual void heatbath_update(std::shared_ptr<SampleState> phi_state,
                               const unsigned int ell) {
    mpi_parallel::cerr
        << "ERROR: heat bath update not implemented for this action "
        << std::endl;
    mpi_exit(EXIT_FAILURE);
  }

  /** @brief Perform local overrelaxation update
   *
   * Update the local entry at position ell of the state using overrelaxation.
   * Usually the index ell corresponds to a dof-index, but depending on
   * the action it can also be the index of a topological entity.
   *
   *  @param[inout] phi_state State to update
   *  @param[in] ell index to update
   */
  virtual void overrelaxation_update(std::shared_ptr<SampleState> phi_state,
                                     const unsigned int ell) {
    mpi_parallel::cerr
        << "ERROR: overrelaxation update not implemented for this action "
        << std::endl;
    mpi_exit(EXIT_FAILURE);
  }

  /** @brief Return heatbath index set */
  const std::vector<unsigned int> &get_heatbath_indexset() const {
    return heatbath_indexset;
  }

  /** @brief Calculate force for HMC integrator for a specific state
   *
   * Calculate \f$P = \frac{\partial S[\phi]}{\partial \phi}\f$ for a specific
   * state and return the resulting force as a state.
   *
   * @param[in] phi_state State \f$\phi\f$ on which to evaluate the force
   * @param[out] p_state Resulting force \f$P\f$ at every point
   *
   */
  virtual void force(const std::shared_ptr<SampleState> phi_state,
                     std::shared_ptr<SampleState> p_state) const = 0;

  /** @brief Initialise state
   *
   * Set initial values of state, those values will be used to start the
   * sampling process
   *
   * @param[out] phi_state State \f$\phi\f$ to be set
   */
  virtual void
  initialise_state(std::shared_ptr<SampleState> phi_state) const = 0;

  /** @brief Copy coarse data points from sample on coarser level
   *
   * @param[in] phi_coarse Coarse sample to copy from
   * @param[in] phi_state Fine state to copy to (sample level as action)
   */
  virtual void copy_from_coarse(const std::shared_ptr<SampleState> phi_coarse,
                                std::shared_ptr<SampleState> phi_state) = 0;

  /** @brief Copy coarse data points from state on finer level
   *
   * @param[in] phi_fine Fine state to copy from
   * @param[in] phi_coarse Coarse state to copy to (same level as action)
   */
  virtual void copy_from_fine(const std::shared_ptr<SampleState> phi_fine,
                              std::shared_ptr<SampleState> phi_state) = 0;

  /** @brief Get coarsening level
   *
   * This will return the coarsening level of the underlying lattice */
  virtual int get_coarsening_level() const = 0;

  /** @brief Action information string
   *
   * return some information on this instance of the action
   */
  virtual std::string info_string() const = 0;

protected:
  /** @brief Renormalisation */
  const RenormalisationType renormalisation;
  /** @brief Set of overrelaxed heatbath indices
   *
   * The OverrelaxedHeatBathSampler will only iterate over these indices,
   * which normally represent dof-indices, but could also correspond
   * to lattice entities if there is more than one dof per entity.
   * If this vector has length zero, all dof-indices will be iterated over.
   */
  std::vector<unsigned int> heatbath_indexset;
};

#endif // ACTION_HH
