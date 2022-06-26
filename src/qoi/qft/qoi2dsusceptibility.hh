#ifndef QOI2DSUSCEPTIBILITY_HH
#define QOI2DSUSCEPTIBILITY_HH QOI2DSUSCEPTIBILITY_HH
#include "action/qft/qftaction.hh"
#include "common/auxilliary.hh"
#include "common/samplestate.hh"
#include "lattice/lattice2d.hh"
#include "mpi/mpi_wrapper.hh"
#include "qoi/quantityofinterest.hh"
#include <cmath>
#include <memory>

/** @file qoi2dsusceptibility.hh
 * @brief Header file for topological susceptibility in the 2D Schwinger model
 */

/** @class QoI2DSusceptibility
 *
 * @brief class for calculating the (dimensionless) topological susceptibility
 * in the 2D Schwinger model
 *
 * The returned topological susceptibility is defined as
 *
 * \f[
 *    V\chi_t = \langle Q^2\rangle
 * \f]
 *
 * where \f$Q\f$ is the topological charge. With this definition \f$V\chi_t\f$
 * approaches the continuum values as \f$\beta\rightarrow\infty\f$ while
 * \f$\beta/P\f$ is kept fixed, with \f$P\f$ denoting the number of plaquettes.
 *
 */

class QoI2DSusceptibility : public QoI {
public:
  /** @brief Create new instance
   *
   * @param[in] lattice_ Lattice
   */
  QoI2DSusceptibility(const std::shared_ptr<Lattice2D> lattice_)
      : lattice(lattice_), four_pi2_inv(0.25 / (M_PI * M_PI)),
        Mt_lat(lattice_->getMt_lat()), Mx_lat(lattice_->getMx_lat()) {}

  /** @brief Destructor */
  virtual ~QoI2DSusceptibility() {}

  /** @brief Evaluate on a state
   *
   * @param[in] phi_state State \f$\phi\f$ on which to evaluate the QoI
   */
  const double virtual evaluate(const std::shared_ptr<SampleState> phi_state);

private:
  /** @brief Temporal extent of lattice */
  const std::shared_ptr<Lattice2D> lattice;
  /** @brief Scaling factor \f$1/(4\pi^2)\f$ */
  const double four_pi2_inv;
  /** @brief Number of time slices */
  const unsigned int Mt_lat;
  /** @brief Number of lattice points in spatial direction */
  const unsigned int Mx_lat;
};

/** @class QoI2DSusceptibilityFactory
 *
 * @brief Factory for constructing the QoI for a particular action
 */
class QoI2DSusceptibilityFactory : public QoIFactory {
public:
  /** @brief Return QoI for a specific  action
   *
   * @param[in] action Action to use
   */
  virtual std::shared_ptr<QoI> get(std::shared_ptr<Action> action) {
    std::shared_ptr<Lattice2D> lattice =
        std::dynamic_pointer_cast<QFTAction>(action)->get_lattice();
    return std::make_shared<QoI2DSusceptibility>(lattice);
  }
};

/** @brief Analytical result for topological susceptibility scaled by the volume
 *
 * Computes the analytical value of the topological susceptibility times the
 * lattice volume, i.e. \f$V\chi_t\f$. If the number of plaquettes is \f$P\f$,
 * then this is given by
 *
 * \f[
 *   V \chi_t(\beta,P) = \langle Q^2 \rangle
 *    = P/\beta \Phi_{\chi_t}(\beta,P)
 * \f]
 *
 * with the function \f$\Phi_{\chi_t}(\beta,P)\f$ define in qoi/auxilliary.hh.
 *
 * For more details see the following two papers:
 *
 *  - Bonati, C. and Rossi, P., 2019. "Topological susceptibility of
 * two-dimensional U (N) gauge theories." Physical Review D, 99(5), p.054503
 * https://doi.org/10.1103/PhysRevD.99.054503
 *
 *  - Kiskis, J., Narayanan, R. and Sigdel, D., 2014. "Correlation between
 * Polyakov loops oriented in two different directions in SU(N) gauge theory on
 * a two-dimensional torus. Physical Review D, 89(8), p.085031.
 *   https://doi.org/10.1103/PhysRevD.89.085031
 *
 *   @param[in] beta Coupling constant \f$beta\f$
 *   @param[in] n_plaq Number of plaquettes \f$P\f$
 */
double quenchedschwinger_chit_analytical(const double beta,
                                         const unsigned int n_plaq);

/** @brief Perturbative approximation to analytical result for topological
 * susceptibility
 *
 * Compute the topological susceptibility up to (and including) corrections of
 * \f$O(1/\beta)\f$. Use this formula if the analytical expression can not be
 * computed reliably due to numerical instabilities which occur if \f$\beta\gg
 * 1\f$.
 *
 *   @param[in] beta Coupling constant \f$beta\f$
 *   @param[in] n_plaq Number of plaquettes \f$P\f$
 */
double quenchedschwinger_chit_perturbative(const double beta,
                                           const unsigned int n_plaq);

/** @brief Analytical result for the variance of the topological susceptibility
 * (scaled by the volume) in the continuum limit
 *
 * Computes the analytical value of the variance of the topological
 * susceptibility times the lattice volume, in the continuum limit i.e.
 * \f$\lim_{\beta\rightarrow \infty}Var[V\chi_t]\f$.
 *
 * If the number of plaquettes is \f$P\f$ and the coupling constant is
 * \f$\beta\f$, then this is given by
 *
 * \f[
 *   V \chi_t(\beta/P) = \left(\hat{\Sigma}_4(\zeta) -
 * \hat{\Sigma}_2(\zeta)^2\right) \f]
 *
 * with \f$\zeta = 4\pi^2\frac{\beta}{P}\f$ and the sums
 *
 *  \f[
 *   \begin{aligned}
 *     \hat{\Sigma}_p(\zeta) &= \frac{\Sigma_p(\zeta)}{\Sigma_0(\zeta)}\\
 *     \Sigma_p(\zeta) &= \sum_{m\in\mathbb{Z}}m^p\exp\left[-\frac{1}{2}\zeta
 * m^2\right] \end{aligned} \f]
 *
 * Note that the limit is taken such that the physical volume is constant, i.e.
 * \f$P/\beta\f$ is fixed.
 *
 *   @param[in] beta Coupling constant \f$beta\f$
 *   @param[in] n_plaq Number of plaquettes \f$P\f$
 */
double
quenchedschwinger_var_chit_continuum_analytical(const double beta,
                                                const unsigned int n_plaq);

#endif // QOI2DSUSCEPTIBILITY_HH
