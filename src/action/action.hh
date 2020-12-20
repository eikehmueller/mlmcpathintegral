#ifndef ACTION_HH
#define ACTION_HH ACTION_HH
#include <memory>
#include <cassert>
#include <random>
#include <vector>
#include <iostream>
#include "common/samplestate.hh"
#include "action/renormalisation.hh"
#include "mpi/mpi_wrapper.hh"

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
    Action(const RenormalisationType renormalisation_) :
        renormalisation(renormalisation_) {}

    /** @brief return size of samples */
    virtual unsigned int sample_size() = 0;

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
    virtual const double evaluate(const std::shared_ptr<SampleState> phi_state) const = 0;

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
    virtual void initialise_state(std::shared_ptr<SampleState> phi_state) const = 0;

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


    /** @brief Check whether action supports number of coarsening steps
     *
     * @param[in] n_level Number of additional coarsening steps (can be zero)
     */
    virtual void check_coarsening_is_permitted(const unsigned int n_level) = 0;

protected:
    /** @brief Renormalisation */
    const RenormalisationType renormalisation;
};

#endif // ACTION_HH
