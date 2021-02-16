#ifndef QOI2DSUSCEPTIBILITY_HH
#define QOI2DSUSCEPTIBILITY_HH QOI2DSUSCEPTIBILITY_HH
#include <memory>
#include <cmath>
#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_integration.h>
#include "common/auxilliary.hh"
#include "common/samplestate.hh"
#include "mpi/mpi_wrapper.hh"
#include "lattice/lattice2d.hh"
#include "action/qft/qftaction.hh"
#include "qoi/quantityofinterest.hh"

/** @file qoi2dsusceptibility.hh
 * @brief Header file for topological susceptibility in the 2D Schwinger model
 */

/** @class QoI2DSusceptibility
 *
 * @brief class for calculating the (dimensionless) topological susceptibility in the 2D Schwinger model
 *
 * The returned topological susceptibility is defined as
 *
 * \f[
 *    \chi_t = \langle Q^2\rangle
 * \f]
 *
 * where \f$Q\f$ is the topological charge. With this definition \f$\chi_t\f$ approaches the
 * continuum values as \f$\beta\rightarrow\infty\f$ while \f$\beta/P\f$ is kept fixed, with
 * \f$P\f$ denoting the number of plaquettes.
 *
 */

class QoI2DSusceptibility : public QoI {
public:
    /** @brief Create new instance
     *
     * @param[in] lattice_ Lattice
     */
    QoI2DSusceptibility(const std::shared_ptr<Lattice2D> lattice_) :
        lattice(lattice_),
        four_pi2_inv(0.25/(M_PI*M_PI)),
        Mt_lat(lattice_->getMt_lat()),
        Mx_lat(lattice_->getMx_lat()) {}

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
        std::shared_ptr<Lattice2D> lattice = std::dynamic_pointer_cast<QFTAction>(action)->get_lattice();
        return std::make_shared<QoI2DSusceptibility>(lattice);
    }
};

/** @brief Analytical result for topological susceptibility scaled by the volume
 *
 * Computes the analytical value of the topological susceptibility times the lattice volume, i.e. \f$V\chi_t\f$.
 * If the number of plaquettes is \f$P\f$, then this is given by
 *
 * \f[
 *   V \chi_t(\beta,P) = \langle Q^2 \rangle
 *    = P\sum_{n=-\infty}^{\infty} w_n(\beta,P)\left[
 *      (P-1) \left(\frac{I'_n(\beta)}{I_n(\beta)}\right)^2 - \frac{I''_n(\beta)}{I_n(\beta)}
 *   \right]
 * \f]
 *  with the weights
 * \f[
 *  w_n(\beta,P) = \frac{(I_n(\beta)^P)}{\sum_{n=-\infty}^{\infty} (I_n(\beta))^P}
 * \f]
 *  and the functions
 *  \f[
 *  \begin{aligned}
 *    I_n(x) &= \frac{1}{2\pi} \int_{-\pi}^{+\pi} e^{i n\phi + x(\cos(\phi)-1)} \; d\phi \\
 *        &= \frac{1}{2\pi} \int_{-\pi}^{+\pi} e^{x(\cos(\phi)-1)} \cos(n\phi) \; d\phi \\
 *    I'_n(x) &= \frac{i}{2\pi} \int_{-\pi}^{+\pi} \frac{\phi}{2\pi }e^{i n\phi + x(\cos(\phi)-1)} \; d\phi \\
 *        &= -\frac{1}{4\pi^2} \int_{-\pi}^{+\pi} \phi e^{x(\cos(\phi)-1)} \sin(n\phi) \; d\phi \\
 *    I''_n(x) &= \frac{i}{2\pi} \int_{-\pi}^{+\pi} \left(\frac{\phi}{2\pi }\right)^2 e^{i n\phi + x(\cos(\phi)-1)} \; d\phi \\
 *        &= -\frac{1}{8\pi^3} \int_{-\pi}^{+\pi} \phi^2 e^{x(\cos(\phi)-1)} \cos(n\phi) \; d\phi
 *  \end{aligned}
 *  \f]
 *  All functions are scaled by a factor \f$e^{-x}\f$ for numerical stability (this factor will cancel out, since we
 *  only every compute ratios of the above functions). Note that \f$I_n(x)\f$ is the (rescaled) modified Bessel
 *  function of the first kind. For more details see the following two papers:
 *
 *  - Bonati, C. and Rossi, P., 2019. "Topological susceptibility of two-dimensional U (N) gauge theories."
 *   Physical Review D, 99(5), p.054503 https://doi.org/10.1103/PhysRevD.99.054503
 *
 *  - Kiskis, J., Narayanan, R. and Sigdel, D., 2014. "Correlation between Polyakov loops oriented in two different
 *   directions in SU(N) gauge theory on a two-dimensional torus. Physical Review D, 89(8), p.085031.
 *   https://doi.org/10.1103/PhysRevD.89.085031
 *
 *   @param[in] beta Coupling constant \f$beta\f$
 *   @param[in] n_plaq Number of plaquettes \f$P\f$
 */
double quenchedschwinger_chit_analytical(const double beta,
                                         const unsigned int n_plaq);

/** @brief Compute functions required in calculation of topological susceptibility
 *
 * Evaluate functions \f$I_n(x)\f$, \f$I'_n(x)\f$ and \f$I''_n(x)\f$
 * required in calculation of \f$\chi_t\f$.
 *
 * @param[in] x Value at which to evaluate functions
 * @param[inout] In Vector which will contain \f$I_n(\beta)\f$
 * @param[inout] dIn Vector which will contain \f$I'_n(\beta)\f$
 * @param[inout] ddIn Vector which will contain \f$I''_n(\beta)\f$
 */
void compute_In(const double x,
                std::vector<double>& In,
                std::vector<double>& dIn,
                std::vector<double>& ddIn);

#endif // QOI2DSUSCEPTIBILITY_HH
