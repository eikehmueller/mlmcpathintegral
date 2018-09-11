#ifndef CLUSTERACTION_HH
#define CLUSTERACTION_HH CLUSTERACTION_HH
#include "path.hh"
#include "action.hh"

/** @file clusteraction.hh
 * @brief Header file for abstract cluster action class
 */

/** @class Action
 *
 * @brief Base class for cluster action
 *
 * Extends the action class by methods which are required to implement the
 * cluster algorithm for one-dimensional quantum problem with periodic
 * boundary conditions. This assumes that the energy can be written as the
 * sum over links
 * \f[
 *   E = \sum_{\ell} E_{\ell}(x^{(\ell)}_-,x^{(\ell)}_+)
 * \f] 
 * where a link \f$\ell=(i,i+1)\f$ connects two neighbouring sites.
 * This, of course, implies that \f$x^{(\ell)}_-=x_i\f$ and
 * \f$x^{(\ell)}_+=x_{i+1}\f$.
 */
class ClusterAction : public Action {
public:
  /** @brief Initialise class
   *
   * Create new instance of class.
   * 
   * @param[in] M_lat_ Number of time slices \f$M\f$
   * @param[in] renormalisation_ Renormalisation type
   * @param[in] T_final_ Final time \f$T\f$
   * @param[in] m0_ Mass of particle \f$m_0\f$
   */
  ClusterAction(const unsigned int M_lat_,
                const double T_final_,
                const RenormalisationType renormalisation_,
                const double m0_)
    : Action(M_lat_, T_final_, renormalisation_, m0_) {}

  /** @brief Change \f$S_{\ell}\f$ in energy used in bonding probabilities
   *
   * The probability to have a bond between sites \f$i\f$ and \f$i+1\f$
   * is given by \f$1-e^{\min(0,-S_{\ell})}\f$
   *
   * @param[in] x_m Value of \f$x^{(\ell)}_- = x_i\f$
   * @param[in] x_p Value of \f$x^{(\ell)}_+ = x_{i+1}\f$
   */
  virtual double S_ell(const double x_m, const double x_p) const = 0;

  /** @brief Set angle for the next step of the cluster algorithm
   */
  virtual void new_angle() const = 0;

  /** @brief Flip the current spin
   * 
   * Return \f$hx\f$
   *
   * @param[in] x value of site \f$x\f$
   */
  virtual double flip(const double x) const = 0;
};

#endif // CLUSTERACTION_HH
