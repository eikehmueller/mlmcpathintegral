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
   * @param[in] T_final_ Final time \f$T\f$
   * @param[in] m0_ Mass of particle \f$m_0\f$
   */
  ClusterAction(const unsigned int M_lat_,
                const double T_final_,
                const double m0_)
    : Action(M_lat_, T_final_, m0_) {}

  /** @brief Energy of link \f$E_{\ell}\f$
   *
   * @param[in] x_m Value of \f$x^{(\ell)}_- = x_i\f$
   * @param[in] x_p Value of \f$x^{(\ell)}_+ = x_{i+1}\f$
   */
  virtual double link_E(const double x_m, const double x_p) const = 0;

  /** @brief Limit on energy of link \f$Q_{\ell}\f$
   *
   * This returns
   * \f[
   *   Q_\ell = \max_{h_-,h_+\in H} E(h_-x^{(\ell)}_-),h_+x^{(\ell)}_+)
   * \f]
   * where \f$H\f$ is the current subgroup. 
   * 
   * @param[in] x_m Value of \f$x^{(\ell)}_- = x_i\f$
   * @param[in] x_p Value of \f$x^{(\ell)}_+ = x_{i+1}\f$
   */
  virtual double link_Q(const double x_m, const double x_p) const = 0;

  /** @brief Set subgroup \f$H\f$ for the next step of the cluster algorithm
   */
  virtual void new_subgroup() const = 0;

  /** @brief Return element \f$h\in H\f$ of currently set subgroup 
   */
  virtual void new_subgroup_element() const = 0;

  /** @brief Apply currently set subgroup \f$h\in H\f$ element to \f$x\f$
   * 
   * Return \f$hx\f$
   *
   * @param[in] x value of site \f$x\f$
   */
  virtual double apply_current_subgroup_element(const double x) const = 0;
};

#endif // CLUSTERACTION_HH
