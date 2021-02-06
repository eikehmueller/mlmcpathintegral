#ifndef QUENCHEDSCHWINGERRENORMALISATION_HH
#define QUENCHEDSCHWINGERRENORMALISATION_HH QUENCHEDSCHWINGERRENORMALISATION_HH
#include "config.h"
#include <memory>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_roots.h>
#include "lattice/lattice2d.hh"
#include "action/renormalisation.hh"
#include "qoi/qft/qoi2dsusceptibility.hh"

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
            betacoarse = betacoarse_nonperturbative();
            break;
        }
        return betacoarse;
    }

private:
    
    /** @brief Non-perturbative value of coarse level coupling
     *
     * The non-perturbative value of the coarse level coupling \f$\beta^{(c)}\f$ is found by
     * matching the topological susceptibility on the fine and coarse level, i.e.
     * \f$\chi_t(\beta^{(c)},P/4)=\chi_t(\beta,P)\f$ where \f$P\f$ is the number of plaquettes
     * on the current level
     */
    double betacoarse_nonperturbative();
    
    /** @brief Type for storing passing function parameters during root finding */
    struct ParamType {
        double beta; // Coupling constant beta
        unsigned int n_plaq; // Number of plaquettes
    };
    
    /** @brief Function used for root finding in non-perturbative matching
     *
     * The non-perturbatively matched coarse level coupling \f$\beta^{(c)}\f$ is
     * found by finding the zeros of the function \f$F(x;\beta,P)= \\chi_t(x\beta,P/4)-\chi_t(\beta,P)\f$
     *
     * @param[in] x Value at which the function is evaluated
     * @param[in] p Parameters (will be cast to ParamType)
     */
    static double f_root (double x, void *p) {
        struct ParamType *params = (struct ParamType *) p;
        double beta = params->beta;
        unsigned int n_plaq = params->n_plaq;
        return  quenchedschwinger_chit_exact(x*beta,n_plaq/4)-quenchedschwinger_chit_exact(beta,n_plaq);
    };
    
    /** @brief Underlying lattice */
    const std::shared_ptr<Lattice2D> lattice;
    /** @brief Coupling constant \f$\beta\f$ */
    const double beta;
};

#endif // QUENCHEDSCHWINGERRENORMALISATION_HH
