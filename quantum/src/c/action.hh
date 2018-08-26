#ifndef ACTION_HH
#define ACTION_HH ACTION_HH
#include <memory>
#include <cassert>
#include <random>
#include <vector>
#include <Eigen/Dense>
#include "path.hh"
#include "sampler.hh"

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
   * @param[in] M_lat_ Number of time slices \f$M\f$
   * @param[in] T_final_ Final time \f$T\f$
   * @param[in] m0_ Mass of particle \f$m_0\f$
   */
  Action(const unsigned int M_lat_,
         const double T_final_,
         const double m0_)
    : M_lat(M_lat_), T_final(T_final_), m0(m0_), a_lat(T_final_/M_lat_) {
    assert(m0>0.0);
    assert(T_final>0.0);
  }

  /** @brief Return number of timeslices \f$M\f$ */
  unsigned int getM_lat() const { return M_lat; }

  /** @return final time \f$T\f$ */
  double getT_final() const { return T_final; }

  /** @brief Return mass \f$m_0\f$ */
  double getm0() const { return m0;}

  /** @brief Return lattice spacing \f$a\f$ */
  double geta_lat() const { return a_lat;}

  /** @brief Evaluate action for a specific path
   * 
   * Calculate \f$S[X]\f$ for a specific path
   *
   * @param[in] x_path Path, has to be an array of length \f$M\f$
   */
  const double virtual evaluate(const std::shared_ptr<Path> x_path) const = 0;

  /** @brief Calculate force for HMC integrator for a specific path
   *
   * Calculate \f$P = \frac{\partial S[X]}{\partial X}\f$ for a specific
   * path and return the resulting force as a path.
   *
   * @param[in] x_path Path \f$X\f$ on which to evaluate the force
   * @param[out] p_path Resulting force \f$P\f$ at every point
   *
   */
  void virtual force(const std::shared_ptr<Path> x_path,
                     std::shared_ptr<Path> p_path) const = 0;

  /** @brief Initialise path 
   *
   * Set initial values of path, those values will be used to start the 
   * sampling process
   *
   * @param[out] x_path Path \f$X\f$ to be set
   */
  void virtual initialise_path(std::shared_ptr<Path> x_path) const = 0;
  
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
  /** @brief Number of time slices */
  const unsigned int M_lat;
  /** @brief Final time */
  const double T_final;
  /** @brief Particle mass */
  const double m0;
  /** @brief lattice spacing */
  const double a_lat;
};

#endif // ACTION_HH
