#ifndef ACTION_HH
#define ACTION_HH ACTION_HH
#include <memory>
#include <cassert>
#include <random>
#include <vector>
#include <iostream>
#include <Eigen/Dense>
#include "lattice/lattice1d.hh"
#include "fields/path.hh"
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
  
   * @param[in] m0_ Mass of particle \f$m_0\f$
   */
  Action(const std::shared_ptr<Lattice1D> lattice_,
         const RenormalisationType renormalisation_,
         const double m0_)
    : lattice(lattice_),
      renormalisation(renormalisation_),
      m0(m0_) {
    assert(m0>0.0);
  }

  /** @brief Return underlying lattice */
  std::shared_ptr<Lattice1D> get_lattice() const { return lattice; }

  /** @brief Return mass \f$m_0\f$ */
  double getm0() const { return m0;}

  /** @brief Cost of one action evaluation */
  virtual double evaluation_cost() const { return lattice->getM_lat(); }

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
  virtual const double evaluate(const std::shared_ptr<Path> x_path) const = 0;

  /** @brief Calculate force for HMC integrator for a specific path
   *
   * Calculate \f$P = \frac{\partial S[X]}{\partial X}\f$ for a specific
   * path and return the resulting force as a path.
   *
   * @param[in] x_path Path \f$X\f$ on which to evaluate the force
   * @param[out] p_path Resulting force \f$P\f$ at every point
   *
   */
  virtual void force(const std::shared_ptr<Path> x_path,
                     std::shared_ptr<Path> p_path) const = 0;

  /** @brief Initialise path 
   *
   * Set initial values of path, those values will be used to start the 
   * sampling process
   *
   * @param[out] x_path Path \f$X\f$ to be set
   */
  virtual void initialise_path(std::shared_ptr<Path> x_path) const = 0;
  
  /** @brief Second derivative \f$W''_{x_-,x_+}(x)\f$ of conditioned
   * action at its minimum.
   *
   * The second derivative (=curvature) of the conditioned action
   * \f$W_{x_-,x_+}(x)=S(\dots,x_-,x,x_+,\dots)\f$ at its minimum. In the
   * special case of a Lagrangian of the form \f$\frac{m_0}{2}\dot{x}^2+V(x)\f$
   * this becomes
   \f[ 
     W_{\overline{x}}(x)=\frac{m_0}{2a}\left((x-x_+)^2+(x-x_-)^2\right)+aV(x)
   \f]
   * where \f$\overline{x}=\frac{x_++x_-}{2}\f$.
   * This quantity is required for sampling of the fine lattice sites.
   *
   * @param[in] x_m Value of \f$x_-\f$
   * @param[in] x_p Value of \f$x_+\f$
   */
  double virtual inline getWcurvature(const double x_m,
                                      const double x_p) const = 0;

  /** @brief Find minimum of conditioned action \f$W_{x_-,x_+}(x)\f$
   *
   * Given \f$x_-\f$ and \f$x_+\f$, find the minimum \f$x_0\f$ 
   * of the conditioned action \f$W_{\overline{x}}(x)\f$
   *
   * @param[in] x_m Value of \f$x_-\f$
   * @param[in] x_p Value of \f$x_+\f$
   */
  double virtual inline getWminimum(const double x_m,
                                    const double x_p) const = 0;
    
protected:
  /** @brief Underlying lattice */
  const std::shared_ptr<Lattice1D> lattice;
  /** @brief Renormalisation */
  const RenormalisationType renormalisation;
  /** @brief Particle mass */
  const double m0;
};

#endif // ACTION_HH
