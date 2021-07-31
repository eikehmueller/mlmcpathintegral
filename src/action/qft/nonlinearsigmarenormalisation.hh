#ifndef NONLINEARSIGMARENORMALISATION_HH
#define NONLINEARSIGMARENORMALISATION_HH NONLINEARSIGMARENORMALISATION_HH
#include "config.h"
#include <memory>
#include "mpi/mpi_wrapper.hh"
#include "lattice/lattice2d.hh"
#include "action/renormalisation.hh"

/** @file nonlinearsigmarenormalisation.hh
 * @brief Header file for renormalisation of nonlinear sigma model
 */

/** @class RenormalisedNonlinearSigmaParameters
 * @brief Renormalised coarse grid parameters for the nonlinear sigma model action
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
     *              (0: none, 1: perturbative [not implemented], 2: nonperturbative
     */
    RenormalisedNonlinearSigmaParameters(const std::shared_ptr<Lattice2D> lattice_,
                                         const double beta_,
                                         const RenormalisationType renormalisation_) :
        RenormalisedParameters(renormalisation_),
        lattice(lattice_),
        beta(beta_) {}

    /** @brief Renormalised coarse level coupling \f$\beta^{(c)}\f$*/
    double beta_coarse() {
        double betacoarse=beta;
        switch (renormalisation) {
            case RenormalisationNone:
                break;
            case RenormalisationPerturbative:
                break;
            case RenormalisationNonperturbative:
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
