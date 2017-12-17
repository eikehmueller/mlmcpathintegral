#ifndef QUARTICOSCILLATORACTION_HH
#define QUARTICOSCILLATORACTION_HH QUARTICOSCILLATORACTION_HH
#include <vector>
#include "path.hh"
#include "action.hh"

/** @file quarticoscillatoraction.hh
 * @brief Header file for quartic oscillator action class
 */


/** @class QuarticOscillatorAction
 *
 * @brief Action class for quartic oscillator 
 *
 * Action class for potential
 * \f$V(x)=\frac{m_0}{2}\mu^2x^2+\frac{\lambda}{4}x^4\f$
 */
class QuarticOscillatorAction : public Action {
public:
  /** @brief Initialise class
   *
   * 
   * @param[in] M_lat_ Number of time slices \f$M\f$
   * @param[in] T_final_ Final time \f$T\f$
   * @param[in] m0_ Mass of particle \f$m_0\f$
   * @param[in] mu2_ Frequency \f$\mu^2\f$
   * @param[in] lambda_ Coefficient of quartic term, \f$\lambda\f$
   */
  QuarticOscillatorAction(const unsigned int M_lat_,
                          const double T_final_,
                          const double m0_,
                          const double mu2_,
                          const double lambda_)
    : Action(M_lat_,T_final_,m0_), mu2(mu2_), lambda(lambda_) {
  }

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
         P_j = \frac{m_0}{a}\left(2X_j-X_{j-1}-X_{j+1}\right) + \frac{1}{2}am_0\mu^2 X_j + \lambda X_j^3
     \f]
   *
   * @param x_path Path \f$X\f$ on which to evaluate the force
   * @param p_path Resulting force \f$P\f$ at every point
   */
  void virtual force(const Path* x_path,
                     Path* p_path) const;

private:
  /** @brief Oscillator frequency */
  const double mu2;
  /** @brief Coefficient of quartic term in potential */
  const double lambda;
};

#endif // QUARTICOSCILLATORACTION_HH
