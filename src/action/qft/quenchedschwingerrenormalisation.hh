#ifndef QUENCHEDSCHWINGERRENORMALISATION_HH
#define QUENCHEDSCHWINGERRENORMALISATION_HH QUENCHEDSCHWINGERRENORMALISATION_HH
#include "config.h"
#include <memory>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_roots.h>
#include "mpi/mpi_wrapper.hh"
#include "lattice/lattice2d.hh"
#include "action/renormalisation.hh"
#include "qoi/qft/qoi2dsusceptibility.hh"

/** @file quenchedschwingerrenormalisation.hh
 * @brief Header file for renormalisation of quenched Schwinger action
 */

/** @class RenormalisedQuenchedSchwingerParameters
 * @brief Renormalised coarse grid parameters for the quenched Schwinger action
 *
 * Calculate the renormalised coarse grid coupling constant \f$\beta\f$ either perturbatively
 * or through non-perturbative matching. Up to \f$\mathcal{O}(\beta^{-1})\f$ we have:
 *
 * \f[
 *    \beta{(c)} = \frac{1}{4}\left(1+\frac{3}{2\beta}\right)\beta
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
     *              (0: none, 1: perturbative [not implemented], 2: nonperturbative
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
        CoarseningType coarsening_type = lattice->get_coarsening_type();
        if (not ( (coarsening_type == CoarsenBoth) or
                  (coarsening_type == CoarsenTemporal) or
                  (coarsening_type == CoarsenSpatial) or
                  (coarsening_type == CoarsenAlternate) ) ) {
                mpi_parallel::cerr << "ERROR: invalid coarsening type in renormalisation" << std::endl;
                mpi_exit(EXIT_FAILURE);
                throw std::runtime_error("...");
            }
        double betacoarse_raw;
        if (coarsening_type == CoarsenBoth) {
            betacoarse_raw = 0.25*beta;
        } else {
            betacoarse_raw = 0.5*beta;
        }
        switch (renormalisation) {
            case RenormalisationNone:
                betacoarse = betacoarse_raw;
                break;
            case RenormalisationPerturbative:
                if (beta > 4.0) {
                    betacoarse = betacoarse_perturbative();
                } else {
                    betacoarse = betacoarse_raw;
                }
                break;
            case RenormalisationNonperturbative:
                if (beta > 4.0) {
                    betacoarse = betacoarse_nonperturbative();
                } else {
                    betacoarse = betacoarse_raw;
                }
                break;
        }
        return betacoarse;
    }

private:
    
    /** @brief Perturbative value of coarse level coupling
     *
     * The perturbative value of the coarse level coupling \f$\beta^{(c)}\f$ is found by
     * matching the topological susceptibility on the fine and coarse level, i.e.
     * \f$V\chi_t(\beta^{(c)},P/4)=\chi_t(\beta,P)\f$ where \f$P\f$ is the number of plaquettes
     * on the current level, only including terms of up to \f$\mathcal{O}(\beta^{-1})\f$
     */
    double betacoarse_perturbative() const {
        double delta;
        double rho;
        if (lattice->get_coarsening_type()==CoarsenBoth) {
            rho = 0.25;
            delta = 1.5;
        } else {
            rho = 0.5;
            delta = 0.5;
        }
        return rho*(1.+delta/beta)*beta;
    }
    
    /** @brief Non-perturbative value of coarse level coupling
     *
     * The non-perturbative value of the coarse level coupling \f$\beta^{(c)}\f$ is found by
     * matching the topological susceptibility on the fine and coarse level, i.e.
     * \f$V\chi_t(\beta^{(c)},P/4)=V\chi_t(\beta,P)\f$ where \f$P\f$ is the number of plaquettes
     * on the current level
     */
    double betacoarse_nonperturbative();
    
    /** @brief Type for storing passing function parameters during root finding */
    struct ParamType {
        double beta; // Coupling constant beta
        unsigned int n_plaq; // Number of plaquettes
        int rho_refine; // Refinement factor
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
        int rho_refine = params->rho_refine;
        return  quenchedschwinger_chit_analytical(x*beta,n_plaq/rho_refine)-quenchedschwinger_chit_analytical(beta,n_plaq);
    };
    
    /** @brief Underlying lattice */
    const std::shared_ptr<Lattice2D> lattice;
    /** @brief Coupling constant \f$\beta\f$ */
    const double beta;
};

#endif // QUENCHEDSCHWINGERRENORMALISATION_HH
