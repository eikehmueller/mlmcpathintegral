#ifndef NONLINEARSIGMAACTION_HH
#define NONLINEARSIGMAACTION_HH NONLINEARSIGMAACTION_HH
#include <memory>
#include <cassert>
#include <random>
#include <vector>
#include <iostream>
#include <algorithm>
#include <iterator>
#include <Eigen/Dense>
#include <Eigen/Geometry>
#include "common/auxilliary.hh"
#include "common/parameters.hh"
#include "common/samplestate.hh"
#include "mpi/mpi_wrapper.hh"
#include "mpi/mpi_random.hh"
#include "distribution/expsin2distribution.hh"
#include "lattice/lattice2d.hh"
#include "action/action.hh"
#include "action/qft/qftaction.hh"

/** @file nonlinearsigmaaction.hh
 * @brief Header file for action of the O(3) non-linear sigma model in two dimensions
 */

/** @class NonlinearSigmaParameters
 *
 * @brief Class for storing parameters of NonlinearSigmaAction
 *
 * This stores the coupling parameter \f$\beta\f$ of the non-linear sigma model
 */
class NonlinearSigmaParameters : public Parameters {
public:
    /** @brief Construct a new instance */
    NonlinearSigmaParameters() :
        Parameters("nonlinearsigma"),
        beta_(1.0),
        renormalisation_(RenormalisationNone){
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
            } else if (renormalisation_str == "nonperturbative") {
                renormalisation_ = RenormalisationNonperturbative;
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

/** @class NonlinearSigmaAction
 *
 * @brief Class for action of two-dimensional O(3) non-linear sigma model
 *
 * Allows calculation with the action
 * \f[
 *      S[\sigma]= - \frac{1}2} \beta \sum_n \sigma_n \cdot (\sigma_{n+\hat{0}+\sigma_{n-\hat{0}+\sigma_{n+\hat{1}+\sigma_{n-\hat{1})
 * \f]
 *
 * with \f$\sigma_n \in\mathbb{R}^3\f$ and \f$|\sigma_n| = 1\f$
 *
 * This is the lattice version of the continuum action
 *
 * \f[
 *      S_{cont}[\sigma] = \frac{1}{2}\beta\int d^2 \sum_{\mu=0,1} \partial_\mu \sigma(x)\cdot \partial_\mu \sigma(x)
 * \f]
 *
 */

class NonlinearSigmaAction : public QFTAction {
public:
    /** @brief Initialise class
     *
     * @param[in] lattice_ Underlying two-dimensional lattice
     * @param[in] lattice_ Underlying fine level two-dimensional lattice
     * @param[in] renormalisation_ Type of renormalisation
     * @param[in] beta_ non-dimensionalised coupling constant \f$\beta=1/(g^2 a_t a_x)\f$
     */
    NonlinearSigmaAction(const std::shared_ptr<Lattice2D> lattice_,
                         const std::shared_ptr<Lattice2D> fine_lattice_,
                         const RenormalisationType renormalisation_,
                         const double beta_)
        : QFTAction(lattice_,fine_lattice_,CoarsenBoth,renormalisation_),
          beta(beta_),
          uniform_dist(0.,2.*M_PI) {
              engine.seed(2481317);
          }

    /** @brief Return coupling constant \f$beta\f$ */
    double getbeta() const {
        return beta;
    }
        
    /** @brief return size of a sample, i.e. the number of links on the lattice */
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
        double beta_coarse = beta;
        /*
        RenormalisedNonlinearSigmaParameters c_param(lattice,
                                                     beta,
                                                     renormalisation);
        double beta_coarse = c_param.beta_coarse();
         */
        std::shared_ptr<Action> new_action;
        std::shared_ptr<Lattice2D> coarse_lattice = lattice->coarse_lattice(2,2,true);
        
        new_action = std::make_shared<NonlinearSigmaAction>(coarse_lattice,
                                                            lattice,
                                                            renormalisation,
                                                            beta_coarse);
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

private:
    /** @brief Extract 3d-vector at lattice vertex (i,j) from sample state and add it to vector
     *
     * @param[in] phi_state State \f$\phi\f$
     * @param[in] i temporal index
     * @param[in] j spatial index
     * @param[out] sigma unit vector
     */
    void add_sigma(const std::shared_ptr<SampleState> phi_state,
                   const int i,
                   const int j,
                   Eigen::Vector3d& sigma) const {
        double theta = phi_state->data[2*lattice->vertex_cart2lin(i,j)];
        double phi = phi_state->data[2*lattice->vertex_cart2lin(i,j)+1];
        sigma[0] += sin(theta)*cos(phi);
        sigma[1] += sin(theta)*sin(phi);
        sigma[2] += cos(theta);
    }

    /** @brief Extract 3d-vector which is given by the sum of the 3d field vectors of
     * the neighbours of site (i,j)
     *
     * @param[in] phi_state State \f$\phi\f$
     * @param[in] i temporal index
     * @param[in] j spatial index
     * @param[out] Delta_n resulting vector
     */
    Eigen::Vector3d delta_neighbours(const std::shared_ptr<SampleState> phi_state,
                                     const int i,
                                     const int j) const {
        Eigen::Vector3d Delta_n(0.,0.,0.);
        add_sigma(phi_state,i+1,j  ,Delta_n);
        add_sigma(phi_state,i-1,j  ,Delta_n);
        add_sigma(phi_state,i  ,j+1,Delta_n);
        add_sigma(phi_state,i  ,j-1,Delta_n);
        return Delta_n;
    }

    

protected:
    /** @brief Dimensionless coupling constant \f$\beta=1/(a_t a_x g^2)\f$*/
    const double beta;
    /** @brief Random number engine */
    typedef mpi_parallel::mt19937_64 Engine;
    /** @brief Type of Mersenne twister engine */
    mutable Engine engine;
    /** @brief distribution for drawing azimuth angle from heat bath */
    mutable std::uniform_real_distribution<double> uniform_dist;
    /** @brief distribution for drawing altitude angle from heat bath */
    const ExpSin2Distribution exp_sin2_dist;
};

#endif // NONLINEARSIGMAACTION_HH
