#ifndef GFFACTION_HH
#define GFFACTION_HH PHI4ACTION_HH
#include <memory>
#include <cassert>
#include <random>
#include <vector>
#include <iostream>
#include <Eigen/Sparse>
#include "common/auxilliary.hh"
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
        mu2_(1.0), 
        renormalisation_(RenormalisationNone) {
          addKey("mu2",Double,Positive);
          addKey("renormalisation",String);
    }

    /** @brief Read parameters from file
     *
     * @param[in] filename Name of file to read
     */
    int readFile(const std::string filename) {

        int readSuccess = Parameters::readFile(filename);
        if (!readSuccess) {
            mu2_ = getContents("mu2")->getDouble();
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

    /** @brief Return squared non-dimensionalised mass \f$\mu^2\f$ */
    double mu2() const {
        return mu2_;
    }

    /** @brief Return renormalisation */
    RenormalisationType renormalisation() const {
        return renormalisation_;
    }
    
private:
    /** @brief Mass parameter \f$\mu^2\f$ */
    double mu2_;
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
 * with the dimensionless squared mass \f$\mu^2\f$.
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
 */

class GFFAction : public QFTAction, public Sampler {
    friend class GFFConditionedFineAction;
public:
    /** @brief Initialise class
     *
     * @param[in] lattice_ Underlying two-dimensional lattice
     * @param[in] lattice_ Underlying fine level two-dimensional lattice
     * @param[in] mu2_ Mass parameter \f$\mu^2\f$
     */
    GFFAction(const std::shared_ptr<Lattice2D> lattice_,
              const std::shared_ptr<Lattice2D> fine_lattice_,
              const double mu2_)
        : QFTAction(lattice_,fine_lattice_,RenormalisationNone),
          mu2(mu2_), sigma(1./sqrt(4.+mu2_)),normal_dist(0.0,1.0) {
              engine.seed(2481317);
              rhs_sample.resize(sample_size());
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
                                                                         2.*mu2);
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
    
    /** @brief Action information string
     *
     * return some information on this instance of the action
     */
    virtual std::string info_string() const;
    
    /** @brief Build (sparse) Cholesky decomposition for exact sampler
     */
    void buildCholesky();
    
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
        
protected:
    /** @brief Squared parameter \f$\mu^2\f$*/
    const double mu2;
    /** @brief Width of Gaussian \f$\sigma = 1\sqrt{4+\mu^2}\f$ */
    const double sigma;
    /** @brief Random number engine */
    mutable mpi_parallel::mt19937_64 engine;
    /** @brief Distribution for heat-bath update */
    mutable std::normal_distribution<double> normal_dist;
    /** @brief Sparse Cholesky matrix L^T for direct sampling */
    mutable Eigen::SparseMatrix<double> choleskyLT;
    /** @brief Vector used for direct sampling*/
    mutable Eigen::VectorXd rhs_sample;
};

#endif // GFFACTION_HH
