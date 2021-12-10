#ifndef GFFACTION_HH
#define GFFACTION_HH PHI4ACTION_HH
#include <memory>
#include <cassert>
#include <random>
#include <vector>
#include <iostream>
#include <Eigen/Sparse>
#include "common/parameters.hh"
#include "common/samplestate.hh"
#include "mpi/mpi_wrapper.hh"
#include "mpi/mpi_random.hh"
#include "lattice/lattice2d.hh"
#include "action/action.hh"
#include "action/qft/qftaction.hh"
#include "sampler/sampler.hh"

/** @file gffaction.hh
 * @brief Header file for action of the 2d Gaussian Free Field (GFF)
 */

/** @class GFFParameters
 *
 * @brief Class for storing parameters of the GFF action
 *
 * This stores the non-dimensionalised squared mass \f$\mu^2=a^2 m^2_{\text{cont}}\f$
 */
class GFFParameters : public Parameters {
public:
    /** @brief Construct a new instance */
    GFFParameters() :
        Parameters("gff"),
        mass_(1.0), 
        renormalisation_(RenormalisationNone) {
          addKey("mass",Double,Positive);
          addKey("renormalisation",String);
    }

    /** @brief Read parameters from file
     *
     * @param[in] filename Name of file to read
     */
    int readFile(const std::string filename) {

        int readSuccess = Parameters::readFile(filename);
        if (!readSuccess) {
            mass_ = getContents("mass")->getDouble();
            std::string renormalisation_str = getContents("renormalisation")->getString();
            if (renormalisation_str == "none") {
                renormalisation_ = RenormalisationNone;
            } else if (renormalisation_str == "perturbative") {
                renormalisation_ = RenormalisationPerturbative;
            } else if (renormalisation_str == "nonperturbative") {
                renormalisation_ = RenormalisationNonperturbative;
            }
        }
        return readSuccess;
    }

    /** @brief Return mass (in units of the inverse physical lattice size) \f$m\f$ */
    double mass() const {
        return mass_;
    }

    /** @brief Return renormalisation */
    RenormalisationType renormalisation() const {
        return renormalisation_;
    }
    
private:
    /** @brief Mass parameter \f$m\f$ */
    double mass_;
    /** @brief Renormalisation */
    RenormalisationType renormalisation_;
};

/** @class GFFAction
 *
 * @brief Class for action of 2d Gaussian Free Field (GFF)
 *
 * Allows calculation with the action
 * \f[
 *      S[\phi] = \sum_n \frac{1}{2}\left[  (\phi_{n+\hat{0}-\phi_n)(\phi_{n+\hat{0}-\phi_n)
 *                                       +  (\phi_{n+\hat{1}-\phi_n)(\phi_{n+\hat{1}-\phi_n)
 *              + \sum_n \frac{1}{2} \mu^2 \phi_n^2
 * \f]
 *
 * with the dimensionless squared mass \f$\mu^2=a^2\mu^2\f$.
 * We assume that all length are measured in units of the lattice
 * size \f$L\f$. Since the correlation length is \f$m^{-1}\f$, this 
 * means that finite volume effects are suppressed provided \f$m\gg1\f$.
 * 
 * To sample directly from the multivariate normal distribution
 * 
 * \f[
 *    \pi(\phi)\sim \exp\left[-\frac{1}{2}\phi^T Q \phi\right]
 * \f]
 *
 * with precision matrix \f$Q\f$, the Cholesky factor \f$L^T\f$ with
 * \f$Q=L L^T\f$ is stored. Samples \f$\phi\f$ from \f$\pi\f$ can
 * then be drawn by drawing a vector \f$\psi\f$ of uncorrelated
 * normal values with mean 0 and variance 1 and solving 
 * \f$L^T\phi=\psi\f$.
 * 
 * On a given level the action can be seen as an approximation of the
 * effective action ontained by integrating out the 'fine' unknowns on
 * the next finer level. This effective action is described by a nine-point
 * stencil and it can be written down as follows:
 * 
 * \f[
 *      S_{eff}[\phi] = \sum_n \frac{1}{2}\left[ \left(4+\mu_f^2 - \frac{4}{4+\mu_f^2}\right) \phi_n^2
 *                                             - \frac{2}{4+\mu_f^2}\left(\phi_{n+\hat{0}}+\phi_{n-\hat{0}}+\phi_{n+\hat{1}}+\phi_{n-\hat{1}}\right) \phi_n
*                                              - \frac{1}{4+\mu_f^2}\left(\phi_{n+\hat{0}+\hat{1}}+\phi_{n+\hat{0}-\hat{1}}+\phi_{n-\hat{0}+\hat{1}}+\phi_{n-\hat{0}-\hat{1}}\right) \phi_n
 *                                         \right]
 *
 * \f]
 *
 * with the non-dimensionalised fine-level mass \f$\mu^f=\frac{1}{2}\mu\f$.
 * 
 * To approximately sample from this action, one can proceed as follows:
 * 
 *   1. Draw a sample from the action \f$S\f$
 *   2. Apply \f$n_{Gibbs}\f$ steps of a Gibbs sampler with the stencil of \f$S_{eff}\f$
 * 
 * As shown in [1], the resulting sample comes from a normal distribution with
 * mean zero and covariance matrix
 * \f[
 *      \hat{\Sigma} = \Sigma_{eff} + G \left(\Sigma - \Sigma_{\eff}\right) G^T
 * \f]
 * 
 * where \f$\Sigma\f$ and \f$\Sigma_{eff}\f$ are the covariance matrices of
 * the actions \f$S\f$ and \f$S_{eff}\f$ respectively and
 * 
 * \f[ 
 *      G = \left((D+L)^{-1}Q_{eff}-I\right)^{n_{Gibbs}}.
 * \f]
 *
 * Here we assumed that the (symmetric) effective precision matrix can be split
 * into a diagonal part \f$D\f$ and a lower triangular part \f$L\f$ as
 * \f$Q_{eff} = D+L+L^T\f$.
 * 
 * [1] Fox, C. and Parker, A., 2017. Accelerated Gibbs sampling of normal
 * distributions using matrix splittings and polynomials. 
 * Bernoulli, 23(4B), pp.3711-3743.
 */

class GFFAction : public QFTAction, public Sampler {
    friend class GFFConditionedFineAction;
public:
    /** @brief Initialise class
     *
     * @param[in] lattice_ Underlying two-dimensional lattice
     * @param[in] lattice_ Underlying fine level two-dimensional lattice
     * @param[in] mass_ Mass parameter \f$m\f$
     * @param[in] n_gibbs_smooth
     */
    GFFAction(const std::shared_ptr<Lattice2D> lattice_,
              const std::shared_ptr<Lattice2D> fine_lattice_,
              const double mass_,
              const int n_gibbs_smooth_=0)
        : QFTAction(lattice_,fine_lattice_,RenormalisationNone),
          mass(mass_), n_gibbs_smooth(n_gibbs_smooth_), normal_dist(0.0,1.0) {
              if (lattice->getMt_lat() != lattice->getMx_lat()) {
                  mpi_parallel::cerr << "ERROR: Lattice has to be squared for GFF action " << std::endl;
                  mpi_exit(EXIT_FAILURE);
              }
              double a_lat;
              if (lattice->is_rotated()) {
                  a_lat = sqrt(2.)/lattice->getMt_lat();
              } else {
                  a_lat = 1./lattice->getMt_lat();
              }
              mu2 = a_lat*a_lat*mass*mass;
              sigma = 1./sqrt(4.+mu2);
              engine.seed(2481317);
              rhs_sample.resize(sample_size());
              buildMatrices();
          }          

    /** @brief Return squared mass parameter \f$\mu^2\f$ */
    double getmu2() const {
        return mu2;
    }
            
    /** @brief return size of a sample, i.e. the number of sites on the lattice */
    virtual unsigned int sample_size() const {
        return lattice->getNvertices();
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
        // Construct coarse lattice
        std::shared_ptr<Lattice2D> coarse_lattice = lattice->get_coarse_lattice();
        // Construct coarse action based on this lattice
        std::shared_ptr<Action> new_action = std::make_shared<GFFAction>(coarse_lattice,
                                                                         lattice,
                                                                         mass,1);
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
    
    /** @brief Global heat bath update with effective action
     * 
     * Carry out a global Gibbs update of the state, traversing the unknowns in
     * lexicographic order.
     * 
     * @param[inout] phi_state Sample state to be updated
     */
    void global_heatbath_update_eff(std::shared_ptr<SampleState> phi_state) const;
    
    /** @brief Build precision matrix based on a given stencil
     * 
     * The stencil is given as a list of entries, which decsribe the
     * coupling on the diagonal (index 0), the nearest neighbours
     * (index 1) and the diagonal nearest neighbours (index 2).
     * For a five-point stencil, the list contains only two entries,
     * for a nine-point stencil it has three entries.
     * 
     * @param[in] stencil List of stencil entries, can have any length
     */
    Eigen::SparseMatrix<double> buildPrecisionMatrix(std::vector<double> stencil);
    
    /** @brief Action information string
     *
     * return some information on this instance of the action
     */
    virtual std::string info_string() const;
    
    /** @brief Build (sparse) Cholesky decomposition for exact sampler
     * as well as covariance matrices used for Gibbs acceleration.
     */
    void buildMatrices();
    
    /** @brief Draw new sample
     * 
     * If the action is Gaussian (\f$\lambda=0\f$ and \f$\mu_0^2<0\f$),
     * draw a direct sample using the Cholesky-decomposition of the 
     * covariance matrix
     */
     virtual void draw(std::shared_ptr<SampleState> phi_state);
     
     /** @brief Set current sample state
      * 
      * This does not do anything since the sample is exact.
      */
     virtual void set_state(std::shared_ptr<SampleState> phi_state) {};
        
     /** @brief print out correlation function */
     void print_correlation_function() const;
    
protected:
    /** @brief Mass parameter \f$m\f$*/
    const double mass;
    /** @brief Non-dimensionalised mass \f$\mu^2 = a^2m^2\f$*/
    double mu2;
    /** @brief Number of Gibbs smoothing steps */
    const int n_gibbs_smooth;
    /** @brief Width of Gaussian for heat-bath update: 
      * \f$\sigma = 1/\sqrt{1+\mu^2}\f$ */
    double sigma;
    /** @brief Random number engine */
    mutable mpi_parallel::mt19937_64 engine;
    /** @brief Distribution for heat-bath update */
    mutable std::normal_distribution<double> normal_dist;
    /** @brief Sparse Cholesky matrix L^T for direct sampling */
    mutable Eigen::SparseMatrix<double> choleskyLT;
    mutable Eigen::SparseMatrix<double> choleskyL;
    /** @brief Precision matrix after Gibbs smoothing */
    mutable Eigen::MatrixXd Q_precision_hat;
    /** @brief Vector used for direct sampling*/
    mutable Eigen::VectorXd rhs_sample;
};

/** @class GFFSamplerFactory
 * 
 * @brief Sampler factory for exact sampling from GFF distribution
 */
class GFFSamplerFactory : public SamplerFactory {
public:
    /** @brief Create new instance
     *
     * @param[in] param_gff GFF parameters
     */
    GFFSamplerFactory() {}

    /** @brief Destructor */
    virtual ~GFFSamplerFactory() {}

    /** @brief Return sampler for a specific  action
     *
     * @param[in] action Action to sample from
     */
    virtual std::shared_ptr<Sampler> get(std::shared_ptr<Action> action) {
        // All we need to do in this case is cast the action to a GFF action
        return std::dynamic_pointer_cast<GFFAction>(action);
    }
};


#endif // GFFACTION_HH
