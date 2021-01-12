#ifndef QUENCHEDSCHWINGERACTION_HH
#define QUENCHEDSCHWINGERACTION_HH QUENCEDSCHWINGERACTION_HH
#include <memory>
#include <cassert>
#include <random>
#include <vector>
#include <iostream>
#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_integration.h>
#include "common/auxilliary.hh"
#include "common/parameters.hh"
#include "common/samplestate.hh"
#include "mpi/mpi_wrapper.hh"
#include "mpi/mpi_random.hh"
#include "distribution/expsin2distribution.hh"
#include "lattice/lattice2d.hh"
#include "action/action.hh"

/** @file quenchedschwingeraction.hh
 * @brief Header file for action of the quenched 2D Schwinger model
 */

/** @class SchwingerParameters
 *
 * @brief Class for storing parameters of Schwinger action
 *
 * This stores the (physical) coupling parameter \f$\beta\f$ of the 2D Schwinger model
 */
class SchwingerParameters : public Parameters {
public:
    /** @brief Construct a new instance */
    SchwingerParameters() :
        Parameters("schwinger"),
        beta_(1.0),
        renormalisation_(RenormalisationNone) {
        addKey("beta",Double,Positive);
        addKey("renormalisation",String);
    }

    /** @brief Read parameters from file
     *
     * @param[in] filename Name of file to read
     */
    int readFile(const std::string filename) {

        int readSuccess = Parameters::readFile(filename);
        if (!readSuccess) {
            beta_ = getContents("beta")->getDouble();
            std::string renormalisation_str = getContents("renormalisation")->getString();
            if (renormalisation_str == "none") {
                renormalisation_ = RenormalisationNone;
            } else if (renormalisation_str == "perturbative") {
                renormalisation_ = RenormalisationPerturbative;
            } else if (renormalisation_str == "exact") {
                renormalisation_ = RenormalisationExact;
            }
        }
        return readSuccess;
    }

    /** @brief Return coupling parameter \f$\beta\f$ */
    double beta() const {
        return beta_;
    }
    /** @brief Return renormalisation */
    RenormalisationType renormalisation() const {
        return renormalisation_;
    }

private:
    /** @brief Coupling parameter \f$\beta\f$ */
    double beta_;
    /** @brief Renormalisation */
    RenormalisationType renormalisation_;
};

/** @class QuenchedSchwingerAction
 *
 * @brief Class for action of quenched 2D Schwinger model
 *
 * Allows calculation with the pure gauge action
 * \f[
 *      S[U]=\beta \sum_P \left[1-Re(U_P) \right]
 * \f]
 * where \f$P = U_{n,0} U_{n+\hat{0},1} U^\dagger_{n+\hat{1},0}U^\dagger_{0,1}\f$
 * is the plaquette at the site (or lattice point) \f$n=(i,j)\f$.
 * \f$beta=\frac{1}{a_t a_x g^2}\f$ is the dimensionless coupling constant.
 *
 * As the link at site \f$n\f$ in direction \f$\mu\f$ is defined as
 * \f[
 *       U_{n,\mu} = e^{i \theta_{n,\mu}}
 * \f]
 * the action can be expressed in terms of the real variables \f$\theta_{n,\mu}\f$ as
 * \f[
 *      S[\theta]=\beta \sum_P \left[1-\cos\left(\theta_{n,0}+\theta_{n+\hat{0},1}-\theta_{n+\hat{1},0}-\theta_{n,1}\right)) \right]
 * \f]
 *
 *  Here if \f$n=(i,j)\f$ then \f$n+\hat{0}=(i+1,j)\f$, \f$n+\hat{1}=(i,j+1)\f$ etc.
 */

class QuenchedSchwingerAction : public Action {
public:
    /** @brief Initialise class
     *
     * @param[in] lattice_ Underlying two-dimensional lattice
     * @param[in] renormalisation_ Type of renormalisation
     * @param[in] beta_ non-dimensionalised coupling constant \f$\beta=1/(g^2 a_t a_x)\f$
     */
    QuenchedSchwingerAction(const std::shared_ptr<Lattice2D> lattice_,
                            const RenormalisationType renormalisation_,
                            const double beta_)
        : Action(renormalisation_),
        lattice(lattice_),
          beta(beta_),
          engine(2481317),
          Mt_lat(lattice->getMt_lat()),
    Mx_lat(lattice->getMx_lat()) { }

    /** @brief Return coupling constant \f$beta\f$ */
    double getbeta() const {
        return beta;
    }
    
    /** @brief return size of a sample, i.e. the number of links on the lattice */
    virtual unsigned int sample_size() const {
        return lattice->getNedges();
    };

    /** @brief Cost of one action evaluation */
    virtual double evaluation_cost() const {
        return sample_size();
    }

    /** @brief Construct coarsened version of action
     *
     * This returns a coarsened version of the action on the next level
     * of the lattice hierarchy.
     */
    virtual std::shared_ptr<Action> coarse_action() {
        std::shared_ptr<Action> new_action;
        new_action = std::make_shared<QuenchedSchwingerAction>(lattice->coarse_lattice(),
                                                               renormalisation,
                                                               beta);
        return new_action;
    };

    /** @brief Evaluate action for a specific state
     *
     * Calculate \f$S[\phi]\f$ for a specific state
     *
     * @param[in] phi_state Sample state \f$\phi\f$
     */
    virtual const double evaluate(const std::shared_ptr<SampleState> phi_state) const;

    /** @brief Draw local value of state from heat bath
     *
     * Update the local entry at position \f$\ell\f$ of the state using a heat bath defined by the neighbouring sites
     *
     *  @param[inout] phi_state State to update
     *  @param[in] ell index of dof to update
     */
    virtual void heatbath_update(std::shared_ptr<SampleState> phi_state, const unsigned int ell);

    /** @brief Perform local overrelaxation update
     *
     * Update the local entry at position \f$\ell\f$ of the state using overrelaxation
     *
     *  @param[inout] phi_state State to update
     *  @param[in] ell index of dof to update
     */
    virtual void overrelaxation_update(std::shared_ptr<SampleState> phi_state, const unsigned int ell);

    /** @brief Calculate force for HMC integrator for a specific state
     *
     * Calculate \f$P = \frac{\partial S[\phi]}{\partial \phi}\f$ for a specific
     * state and return the resulting force as a state.
     *
     * @param[in] phi_state State \f$\phi\f$ on which to evaluate the force
     * @param[out] p_state Resulting force \f$P\f$ at every point
     *
     */
    virtual void force(const std::shared_ptr<SampleState> phi_state,
                       std::shared_ptr<SampleState> p_state) const;

    /** @brief Copy coarse data points from sample on coarser level
     *
     * @param[in] phi_coarse Coarse sample to copy from
     * @param[in] phi_state Fine state to copy to (sample level as action)
     */
    virtual void copy_from_coarse(const std::shared_ptr<SampleState> phi_coarse,
                                  std::shared_ptr<SampleState> phi_state);

    /** @brief Copy coarse data points from state on finer level
     *
     * @param[in] phi_fine Fine state to copy from
     * @param[in] phi_coarse Coarse state to copy to (same level as action)
     */
    virtual void copy_from_fine(const std::shared_ptr<SampleState> phi_fine,
                                std::shared_ptr<SampleState> phi_state);

    /** @brief Initialise state
     *
     * Set initial values of state, those values will be used to start the
     * sampling process
     *
     * @param[out] phi_state State \f$\phi\f$ to be set
     */
    virtual void initialise_state(std::shared_ptr<SampleState> phi_state) const;
    
    /** @brief Get coarsening level
     *
     * This will return the coarsening level of the underlying lattice */
    virtual int get_coarsening_level() const {
        return lattice->get_coarsening_level();
    }

    /** @brief Check whether action supports number of coarsening steps
     *
     * @param[in] n_level Number of additional coarsening steps (can be zero)
     */
    virtual void check_coarsening_is_permitted(const unsigned int n_level);
    
    /** @brief Exact result for topological susceptibility
     *
     * Computes the exact value of the topological susceptibility. If the number of plaquettes is \f$P\f$,
     * then this is given by
     *
     * \f[
     *   \chi_t(\beta,P) = \frac{\langle Q^2 \rangle}{P}
     *   = \sum_{n=-\infty}^{\infty} w_n(\beta,P)\left[
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
     *    I_n(x) &= \frac{1}{2\pi} \int_{-\pi}^{+\pi} e^{i n\phi + x\cos(\phi)} \; d\phi \\
     *        &= \frac{1}{2\pi} \int_{-\pi}^{+\pi} e^{x\cos(\phi)} \cos(n\phi) \; d\phi \\
     *    I'_n(x) &= \frac{i}{2\pi} \int_{-\pi}^{+\pi} \frac{\phi}{2\pi }e^{i n\phi + x\cos(\phi)} \; d\phi \\
     *        &= -\frac{1}{4\pi^2} \int_{-\pi}^{+\pi} \phi e^{x\cos(\phi)} \sin(n\phi) \; d\phi \\
     *    I''_n(x) &= \frac{i}{2\pi} \int_{-\pi}^{+\pi} \left(\frac{\phi}{2\pi }\right)^2 e^{i n\phi + x\cos(\phi)} \; d\phi \\
     *        &= -\frac{1}{8\pi^3} \int_{-\pi}^{+\pi} \phi^2 e^{x\cos(\phi)} \cos(n\phi) \; d\phi
     *  \end{aligned}
     *  \f]
     *  Note that \f$I_n(x)\f$ is the modified Bessel function of the first kind. For more details see the following two papers:
     *
     *  - Bonati, C. and Rossi, P., 2019. "Topological susceptibility of two-dimensional U (N) gauge theories."
     *   Physical Review D, 99(5), p.054503 https://doi.org/10.1103/PhysRevD.99.054503
     *
     *  - Kiskis, J., Narayanan, R. and Sigdel, D., 2014. "Correlation between Polyakov loops oriented in two different
     *   directions in SU(N) gauge theory on a two-dimensional torus. Physical Review D, 89(8), p.085031.
     *   https://doi.org/10.1103/PhysRevD.89.085031
     */
    double chit_exact() const;

private:
        
    /** @brief Compute width-parameter of exponential distribution
     *
     * Given \f$\theta_+,\theta_-\f$ compute \f$\sigma\f$ such that
     * \f[
     *   \cos(\theta-\theta_+) + \cos(\theta-\theta_-) = -\sigma \sin^2\left(\frac{\theta-\theta_0}{2}\right) + const.
     * \f]
     *
     * Explicitly this is given as
     * \f[
     *   \sigma=4|\cos\left(\frac{\theta_+-\theta_-}{2}\right)|
     * \f]
     *
     * NOTE THAT THE GIVEN FORMULAE ARE ONLY CORRECT IF \f$\theta_+,\theta_-\in[-\pi,+\pi]\f$
     *
     * @param[in] theta_p Angle \f$\theta_+\in[-\pi,+\pi]\f$
     * @param[in] theta_m Angle \f$\theta_-\in[-\pi,+\pi]\f$
     */
    inline double sigma_expsin2(const double theta_p, const double theta_m) const {
        return 4.*fabs(cos(0.5*(theta_p-theta_m)));
    }

    /** @brief Compute shift-parameter of exponential distribution
     *
     * Given \f$\theta_+,\theta_-\f$ compute \f$\theta_0\f$ such that
     * \f[
     *   \cos(\theta-\theta_+) + \cos(\theta-\theta_-) = -\sigma \sin^2\left(\frac{\theta-\theta_0}{2}\right) + const.
     * \f]
     *
     * Explicitly this is given as
     * \f[
     *    \theta_0=\begin{cases}
     *    \frac{1}{2} (\theta_+ + \theta_-) & \text{for $\theta_+-\theta_- \in [-\pi,+\pi]$}\\
     *    \frac{1}{2} (\theta_+ + \theta_-)+\pi \mod [-\pi,+pi) & \text{otherwise$}
     *    \end{cases}
     * \f]
     *
     * NOTE THAT THE GIVEN FORMULAE ARE ONLY CORRECT IF \f$\theta_+,\theta_-\in[-\pi,+\pi]\f$
     *
     * @param[in] theta_p Angle \f$\theta_+\in[-\pi,+\pi]\f$
     * @param[in] theta_m Angle \f$\theta_-\in[-\pi,+\pi]\f$
     */
    inline double theta0_expsin2(const double theta_p, const double theta_m) const {
        return mod_2pi(0.5*(theta_p+theta_m)+(fabs((theta_p-theta_m))>M_PI)*M_PI);
    }
    
    /** @brief Compute staple angles \f$\theta_{n,\mu}^{(+)}\f$ and \f$\theta_{n,\mu}^{(-)}\f$
     *
     * Let \f$\nu=1-mu\f$. Then the staple angles are defined as
     * \f[
     *    \theta_{n,\mu}^{(+)} = \theta_{n,\nu} + \theta_{n+\hat{\nu},\mu} - \theta_{n+\hat{\mu},\nu}
     * \f]
     *  and
     * \f[
     *    \theta_{n,\mu}^{(-)} = \theta_{n-\hat{\nu}+\hat{\mu},\nu} + \theta_{n-\hat{\nu},\mu} - \theta_{n-\hat{\nu},\nu}
     * \f]
     * Note that with those angles the action can be written as a sum over all links
     * \f[
     *    S[U] = \frac{\beta}{4}\sum_{n} \sum_{\mu=0,1} \left[2 - \cos(\theta_{n,\mu}-\theta_{n,\mu}^{(+)})
     *                                           - \cos(\theta_{n,\mu}-\theta_{n,\mu}^{(-)})\right]
     * \f]
     *
     *  @param[in] phi_state State vector, containing the angles \f$\theta_{n,\mu}\f$
     *  @param[in] i Temporal index of coordinate \f$n\f$
     *  @param[in] j Temporal index of coordinate \f$n\f$
     *  @param[in] mu Direction \f$\mu\f$
     *  @param[out] theta_p Resulting value of \f$\theta_{n,\mu}^{(+)}\f$
     *  @param[out] theta_m Resulting value of \f$\theta_{n,\mu}^{(-)}\f$
     */
    void compute_staple_angles(std::shared_ptr<SampleState> phi_state,
                               const int i,
                               const int j,
                               const int mu,
                               double& theta_p,
                               double& theta_m);

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
                    std::vector<double>& ddIn) const;

protected:
    /** @brief Underlying lattice */
    const std::shared_ptr<Lattice2D> lattice;
    /** @brief Dimensionless coupling constant \f$\beta=1/(a_t a_x g^2)\f$*/
    const double beta;
    /** @brief Random number engine */
    typedef mpi_parallel::mt19937_64 Engine;
    /** @brief Type of Mersenne twister engine */
    mutable Engine engine;
    /** @brief Number of time slices */
    const unsigned int Mt_lat;
    /** @brief Number of points in spatial direction */
    const unsigned int Mx_lat;
    /** @brief distribution for drawing from heat bath */
    const ExpSin2Distribution exp_sin2_dist;

};

#endif // QUENCHEDSCHWINGERACTION_HH
