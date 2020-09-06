#ifndef ACTION_HH
#define ACTION_HH ACTION_HH
#include <memory>
#include <cassert>
#include <random>
#include <vector>
#include <iostream>
#include "common/samplestate.hh"
#include "lattice/lattice1d.hh"
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
   * 
   * @param[in] lattice_ Underlying lattice
   * @param[in] renormalisation_ Type of renormalisation
    */
  Action(const std::shared_ptr<Lattice1D> lattice_,
         const RenormalisationType renormalisation_)
    : lattice(lattice_),
      renormalisation(renormalisation_) {}

  /** @brief return size of samples */
  virtual unsigned int sample_size() = 0;

  /** @brief Return underlying lattice */
  std::shared_ptr<Lattice1D> get_lattice() const { return lattice; }

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
     
  /** @brief Evaluate action for a specific path
   * 
   * Calculate \f$S[X]\f$ for a specific path
   *
   * @param[in] x_path Path, has to be an array of length \f$M\f$
   */
  virtual const double evaluate(const std::shared_ptr<SampleState> x_path) const = 0;

  /** @brief Calculate force for HMC integrator for a specific path
   *
   * Calculate \f$P = \frac{\partial S[X]}{\partial X}\f$ for a specific
   * path and return the resulting force as a path.
   *
   * @param[in] x_path Path \f$X\f$ on which to evaluate the force
   * @param[out] p_path Resulting force \f$P\f$ at every point
   *
   */
  virtual void force(const std::shared_ptr<SampleState> x_path,
                     std::shared_ptr<SampleState> p_path) const = 0;

  /** @brief Initialise path 
   *
   * Set initial values of path, those values will be used to start the 
   * sampling process
   *
   * @param[out] x_path Path \f$X\f$ to be set
   */
  virtual void initialise_path(std::shared_ptr<SampleState> x_path) const = 0;
      
  /** @brief Copy coarse data points from sample on coarser level
   *
   * @param[in] x_coarse Coarse sample to copy from
   * @param[in] x_path Fine path to copy to (sample level as action)
   */
  virtual void copy_from_coarse(const std::shared_ptr<SampleState> x_coarse,
                                std::shared_ptr<SampleState> x_path) = 0;
  
  /** @brief Copy coarse data points from path on finer level
   *
   * @param[in] x_fine Fine path to copy from
   * @param[in] x_coarse Coarse sample to copy to (same level as action)
   */
  virtual void copy_from_fine(const std::shared_ptr<SampleState> x_fine,
                              std::shared_ptr<SampleState> x_path) = 0;

  /** @brief Get coarsening level
   * 
   * This will return the coarsening level of the underlying lattice */
   int get_coarsening_level() { return lattice->get_coarsening_level(); }
      
      
  /** @brief Check whether action supports number of coarsening steps 
   * 
   * @param[in] n_level Number of additional coarsening steps (can be zero)
   */
   virtual void check_coarsening_is_permitted(const unsigned int n_level) = 0;
   
protected:
  /** @brief Underlying lattice */
  const std::shared_ptr<Lattice1D> lattice;
  /** @brief Renormalisation */
  const RenormalisationType renormalisation;
};

#endif // ACTION_HH
