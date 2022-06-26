#ifndef NONLINEARSIGMARENORMALISATION_HH
#define NONLINEARSIGMARENORMALISATION_HH NONLINEARSIGMARENORMALISATION_HH
#include "action/renormalisation.hh"
#include "config.h"
#include "lattice/lattice2d.hh"
#include "mpi/mpi_wrapper.hh"
#include <memory>

/** @file nonlinearsigmarenormalisation.hh
 * @brief Header file for renormalisation of nonlinear sigma model
 */

/** @class RenormalisedNonlinearSigmaParameters
 * @brief Renormalised coarse grid parameters for the nonlinear sigma model
 * action
 *
 * Calculate the renormalised coarse grid coupling constant \f$\beta\f$.
 */
class RenormalisedNonlinearSigmaParameters : public RenormalisedParameters {
public:
  /** @brief Create new instance
   *
   * @param[in] lattice_ Underlying lattice
   * @param[in] beta_ Coupling constant \f$\beta\f$
   * @param[in] renormalisation_ Type of renormalisation to use
   *              (0: none, 1: perturbative [not implemented], 2:
   * nonperturbative
   */
  RenormalisedNonlinearSigmaParameters(
      const std::shared_ptr<Lattice2D> lattice_, const double beta_,
      const RenormalisationType renormalisation_)
      : RenormalisedParameters(renormalisation_), lattice(lattice_),
        beta(beta_) {}

  /** @brief Renormalised coarse level coupling \f$\beta^{(c)}\f$
   *
   * So far, only perturbative renormalisation has been implemented.
   *
   * If we write the Lagrangian as 1/(2*g^2)*|d/d_mu n|^2, the
   * \f$\beta\f$-function is given by [Peskin & Schroeder section 13.3, Eq.
   * (13.86)]:
   *
   * \beta(g) = -g^3/(4\pi) + O(g^5)
   *
   * Since
   *
   * \f[
   *      \beta = -\frac{dg}{d\log(a M)}
   * \f]
   *
   * (where \f$M\f$ is a reference scale that cancels out).
   *
   * we find that
   *
   * 1/g(2a)^2 = 1/g(a)^2 - \frac{\log(2)}{2\pi} + O(g)
   *
   */
  double beta_coarse() {
    double betacoarse;
    switch (renormalisation) {
    case RenormalisationNone:
      betacoarse = beta;
      break;
    case RenormalisationPerturbative:
      betacoarse = beta - 0.5 * log(2.) / (2. * M_PI);
      break;
    case RenormalisationNonperturbative:
      mpi_parallel::cerr << "ERROR: non-perturbative renormalisation not "
                            "implemented for non-linear sigma model."
                         << std::endl;
      mpi_exit(EXIT_FAILURE);
      throw std::runtime_error("...");
      break;
    }
    return betacoarse;
  }

  /** @brief Underlying lattice */
  const std::shared_ptr<Lattice2D> lattice;
  /** @brief Coupling constant \f$\beta\f$ */
  const double beta;
};

#endif // NONLINEARSIGMARENORMALISATION_HH
