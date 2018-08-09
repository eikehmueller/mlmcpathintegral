#ifndef ROTORACTION_HH
#define ROTORACTION_HH ROTORACTION_HH
#include <math.h>
#include <random>
#include "path.hh"
#include "action.hh"

/** @file rotoraction.hh
 * @brief Header file for quantum mechanical rotor action class
 */


/** @class HarmonicOscillatorAction
 *
 * @brief Action class for the quantum mechanical rotor descibed in
 *        <a href="https://arxiv.org/abs/1503.05088">arXiv/1503.05088</a>  
 *
 * Action class for free particle moving on a circle of radius \f$R\f$ with 
 * moment of inertia (angular mass) \f$I=m_0R^2\f$. The discrete action is
 * given by
 * \f[
 *   S = \frac{I}{a}\sum_{i=0}^{M_{lat}-1} \left(1-\cos(x_i-x_{i-1})\right)
 * \f]
 * with periodic boundary conditions \f$x_{M_{lat}}=x_0\f$
 */
class RotorAction : public Action {
public:
  /** @brief Initialise class
   *
   * 
   * @param[in] M_lat_ Number of time slices \f$M\f$
   * @param[in] T_final_ Final time \f$T\f$
   * @param[in] m0_ Moment of inertia (angular mass) \f$I\f$
   */
  RotorAction(const unsigned int M_lat_,
              const double T_final_,
              const double m0_)
    : Action(M_lat_,T_final_,m0_) {}

  /** @brief Tidy up
   * 
   * Delete any temporary arrays
   */
  ~RotorAction() {}
  
  /** @brief Evaluate action for a specific path
   * 
   * Calculate \f$S[X]\f$ for a specific path
   *
   * @param[in] x_path path \f$X\f$, has to be am array of length \f$M\f$
   */
  const double virtual evaluate(const Path* x_path) const;

  /** @brief Calculate force for HMC integrator for a specific path
   *
   * Calculate \f$P = \frac{\partial S[X]}{\partial X}\f$ for a specific
   * path and return the resulting force as a path \f$P\f$.
   *
   * Note that for this action we have
     \f[
         P_j = \frac{I}{a}\left(\sin(x_j-x_{j-1})+\sin(x_j-x_{j+1})\right)
     \f]
   *
   * @param x_path Path \f$X\f$ on which to evaluate the force
   * @param p_path Resulting force \f$P\f$ at every point
   */
  void virtual force(const Path* x_path,
                     Path* p_path) const;

  /** @brief Initialise path 
   *
   * Set initial values of path to random values in the range [-pi,pi)
   *
   * @param[out] x_path Path \f$X\f$ to be set
   */
  void virtual initialise_path(Path* x_path) const;
  
  /** @brief Second derivative \f$W''_{x_-,x_+}(x)\f$ of conditioned action
   *
   * For the quantum mechanical rotor the curvature of the modified
   * action (see Action::getWcurvature()) is 
   \f[
     W''_{x_-,x_+} = \frac{I}{a}\left|\cos\left(\frac{x_+-x_-}{2}\right)\right|
   \f]
   *
   * @param[in] x_m Value of \f$x_-\f$
   * @param[in] x_p Value of \f$x_+\f$
   */
  double virtual inline getWcurvature(const double x_m,
                                      const double x_p) const {
    return 2.0*m0/a_lat*fabs(cos(0.5*(x_p-x_m)));
  }

  /** @brief Find minimum of conditioned action \f$W_{x_-,x_+}(x)\f$
   *
   * For the quantum mechanical rotor the minimum of the modified
   * action (see Action::getWminimum()) can be found at
   \f[
      \tan(x_0) = \frac{\sin(x_-)+\sin(x_+)}{\cos(x_-)+\cos(x_+)}
   \f]
   * 
   * @param[in] x_m Value of \f$x_-\f$
   * @param[in] x_p Value of \f$x_+\f$
   */
  double virtual inline getWminimum(const double x_m,
                                    const double x_p) const {
    return atan2(sin(x_p)+sin(x_m),cos(x_p)+cos(x_m));
  }

};

#endif // ROTORACTION_HH
