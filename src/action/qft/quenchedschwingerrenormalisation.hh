#ifndef QUENCHEDSCHWINGERRENORMALISATION_HH
#define QUENCHEDSCHWINGERRENORMALISATION_HH QUENCHEDSCHWINGERRENORMALISATION_HH
#include "config.h"
#include <memory>
#include "lattice/lattice2d.hh"
#include "action/renormalisation.hh"

/** @file quenchedschwingerrenormalisation.hh
 * @brief Header file for renormalisation of quenched Schwinger action
 */

/** @class RenormalisedQuenchedSchwingerParameters
 * @brief Renormalised coarse grid parameters for the quenched Schwinger action
 *
 * Calculate the renormalised coarse grid coupling constant \f$\beta\f$
 *
 * \f[
 *    \beta_0^{(c)} = \frac{1}{4}\beta
 * \f]
 *
 */
class RenormalisedQuenchedSchwingerParameters : public RenormalisedParameters {
public:
    /** @brief Create new instance
     *
     * @param[in] lattice_ Underlying lattice
     * @param[in] beta_ Coupling constant \f$\beta\f$
     * @param[in] renormalisation_ Type of renormalisation to use
     *              (0: none, 1: perturbative, 2: exact [not implemented])
     */
    RenormalisedQuenchedSchwingerParameters(const std::shared_ptr<Lattice2D> lattice_,
                                            const double beta_,
                                            const RenormalisationType renormalisation_) :
        RenormalisedParameters(renormalisation_),
        lattice(lattice_),
        beta(beta_) {}

    /** @brief Renormalised coarse level mass \f$\beta^{(c)}\f$*/
    double beta_coarse() {
        double betacoarse;
        switch (renormalisation) {
        case RenormalisationNone:
            betacoarse = 0.25*beta;
            break;
        case RenormalisationPerturbative:
            betacoarse = 1.0;
            mpi_parallel::cerr << "ERROR: perturbative renormalisation not implemented for quenched Schwinger action " << std::endl;
            mpi_exit(EXIT_FAILURE);
            break;
        case RenormalisationExact:
            betacoarse = 1.0;
            mpi_parallel::cerr << "ERROR: exact renormalisation not implemented for quenched Schwinger action " << std::endl;
            mpi_exit(EXIT_FAILURE);
            break;
        }
        return betacoarse;
    }

private:
    /** @brief Underlying lattice */
    const std::shared_ptr<Lattice2D> lattice;
    /** @brief Coupling constant \f$\beta\f$ */
    const double beta;
};

#endif // QUENCHEDSCHWINGERRENORMALISATION_HH
