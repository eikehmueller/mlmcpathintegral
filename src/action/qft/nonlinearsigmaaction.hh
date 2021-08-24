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
#include "action/qft/nonlinearsigmarenormalisation.hh"

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
 *      S[\sigma]= - \frac{1}2} \beta \sum_n \sigma_n \cdot \Delta_n
 * \f]
 *
 * where \f$\Delta_n = \sigma_{n+\hat{0}+\sigma_{n-\hat{0}+\sigma_{n+\hat{1}+\sigma_{n-\hat{1}\f$
 * with \f$\sigma_n \in\mathbb{R}^3\f$ and \f$|\sigma_n| = 1\f$
 *
 * This is the lattice version of the continuum action
 *
 * \f[
 *      S_{cont}[\sigma] = \frac{1}{2}\beta\int d^2 \sum_{\mu=0,1} \partial_\mu \sigma(x)\cdot \partial_\mu \sigma(x)
 * \f]
 *
 * The action can be either formulated on a normal, unrotated lattice or on a rotated
 * lattice. In both cases the same Cartesion lattice object is used, but in the lattice
 * case only the points marked in the following diagram are used:
 * 
 * 
 *   X--+--X--+--X--+--X
 *   !  !  !  !  !  !  !
 *   +--B--+--B--+--X--+
 *   !  !  !  !  !  !  !
 *   X--+--A--+--X--+--X
 *   !  !  !  !  !  !  !
 *   +--B--+--B--+--X--+
 *   !  !  !  !  !  !  !
 *   X--+--X--+--X--+--X
 * 
 * Note that for a rotated lattice the nearest neighbours of a point are defined as the 
 * 'diagonal' neighbours. For example, in the above diagram the point 'A' has the four
 * nearest neighbours marked with the letter 'B'.
 * 
 */

class NonlinearSigmaAction : public QFTAction {
    // The conditioned fine action class needs access to some private class 
    // methods to avoid code duplication
    friend class NonlinearSigmaConditionedFineAction;
public:
    /** @brief Initialise class
     *
     * @param[in] lattice_ Underlying two-dimensional lattice
     * @param[in] lattice_ Underlying fine level two-dimensional lattice
     * @param[in] renormalisation_ Type of renormalisation
     * @param[in] beta_ non-dimensionalised coupling constant \f$\beta=1/(g^2 a_t a_x)\f$
     * @param[in] rotated_ Is the action formulated on a rotated lattice?
     */
    NonlinearSigmaAction(const std::shared_ptr<Lattice2D> lattice_,
                         const std::shared_ptr<Lattice2D> fine_lattice_,
                         const RenormalisationType renormalisation_,
                         const double beta_,
                         const double rotated_=false)
        : QFTAction(lattice_,fine_lattice_,CoarsenBoth,renormalisation_),
          beta(beta_),
          rotated(rotated_),
          uniform_dist(0.,2.*M_PI) {
              engine.seed(2481317);
          }

    /** @brief Return coupling constant \f$beta\f$ */
    double getbeta() const {
        return beta;
    }
    
    /** @brief return true if the action is formulated on a rotated lattice */
    bool is_rotated() const {
        return rotated;
    }
        
    /** @brief return size of a sample, i.e. the number of vertices on the lattice */
    virtual unsigned int sample_size() const {
        if (rotated) {
            return lattice->getNvertices();
        } else {
            return 2*lattice->getNvertices();
        }
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
        RenormalisedNonlinearSigmaParameters c_param(lattice,
                                                     beta,
                                                     renormalisation);
        double beta_coarse = c_param.beta_coarse();
        std::shared_ptr<Action> new_action;
        std::shared_ptr<Lattice2D> coarse_lattice;
        if (rotated) {
            coarse_lattice = lattice->coarse_lattice(2,2,true);
        } else {
            coarse_lattice = lattice;
        }
        
        new_action = std::make_shared<NonlinearSigmaAction>(coarse_lattice,
                                                            lattice,
                                                            renormalisation,
                                                            beta_coarse,
                                                            not rotated);
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
     * Update the local entry at position \f$\ell\f$ of the state using a heat bath
     * defined by the neighbouring sites.
     *
     *  @param[inout] phi_state State to update
     *  @param[in] ell index of dof to update
     */
    virtual void heatbath_update(std::shared_ptr<SampleState> phi_state, const unsigned int ell);

    /** @brief Perform local overrelaxation update
     *
     * Update the local entry at position \f$\ell\f$ of the state using overrelaxation.
     *
     * For this, draw an angle \f$\tau_n\f$ uniformly from the interval \f$[0,2\pi[\f$ and
     * rotate \f$\sigma_n\f$ around \f$\Delta_n\f$ by this angle.
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
     * Since the action can be written as \f$S[\sigma]= - \frac{1}2} \beta \sum_n \sigma_n \cdot \Delta_n\f$,
     * the force at point \f$n\f$is given by
     * \f[
     * \begin{aligned}
     *   \frac{\partial S}/{\partial\theta} &= -\beta\left\[ ((\Delta_n)_0 \cos(\phi_n)+(\Delta_n)_1\sin(\phi_n)) \cos(\theta_n)
     *                                -(\Delta_n)_2\sin(\theta_n) \right] \\
     *   \frac{\partial S}/{\partial\phi} &= -\beta\left[-(\Delta_n)_0\sin(\phi_n)+(\Delta_n)_1\cos(\phi_n)\right]
     * \end{aligned}
     * \f]
     *
     *
     * @param[in] phi_state State \f$\phi\f$ on which to evaluate the force
     * @param[out] p_state Resulting force \f$P\f$ at every point
     *
     */
    virtual void force(const std::shared_ptr<SampleState> phi_state,
                       std::shared_ptr<SampleState> p_state) const;

    /** @brief Get values of \f$\theta,\phi\f$ stored in dof-vector for specific (i,j)
     * 
     * @param[in] phi_state State vector to extract dofs from
     * @param[in] i Temporal index
     * @param[in] j Spatial index
     * @param[out] theta Value of angle theta
     * @param[out] phi Value of angle phi
     */
    virtual void get_dofs(const std::shared_ptr<SampleState> phi_state,
                          const int i,
                          const int j,
                          double& theta,
                          double& phi) const {
        // Linear index
        unsigned int ell;
        if (rotated) {
            ell = lattice->diag_vertex_cart2lin(i,j);
        } else {
            ell = lattice->vertex_cart2lin(i,j);
        }
        theta = phi_state->data[2*ell];
        phi = phi_state->data[2*ell+1];
    }

    /** @brief Set values of \f$\theta,\phi\f$ stored in dof-vector for specific (i,j)
     * 
     * @param[in] phi_state State vector
     * @param[in] i Temporal index
     * @param[in] j Spatial index
     * @param[in] theta Value of angle theta
     * @param[in] phi Value of angle phi
     */
    virtual void set_dofs(std::shared_ptr<SampleState> phi_state,
                          const int i,
                          const int j,
                          const double theta,
                          const double phi) const {
        // Linear index
        unsigned int ell;
        if (rotated) {
            ell = lattice->diag_vertex_cart2lin(i,j);
        } else {
            ell = lattice->vertex_cart2lin(i,j);
        }
        phi_state->data[2*ell] = theta;
        phi_state->data[2*ell+1] = phi;
    }

    /** @brief Copy coarse data points from sample on coarser level
     *
     * There are two cases, depending on whether the action is formulated on
     * a rotated lattice or not:
     *
     * === Case 1 === we are operating on an un-rotated lattice
     * 
     *     current lattice              coarse lattice
     *   C--F--C--F--C--F--C          C--+--C--+--C--+--C
     *   !  !  !  !  !  !  !          !  !  !  !  !  !  !
     *   F--C--F--C--F--C--F          +--C--+--C--+--C--+
     *   !  !  !  !  !  !  !          !  !  !  !  !  !  !
     *   C--F--C--F--C--F--C          C--+--C--+--C--+--C
     *   !  !  !  !  !  !  !          !  !  !  !  !  !  !
     *   F--C--F--C--F--C--F          +--C--+--C--+--C--+
     *   !  !  !  !  !  !  !          !  !  !  !  !  !  !
     *   C--F--C--F--C--F--C          C--+--C--+--C--+--C
     * 
     *  => Copy the C-points from the next coarser lattice (right)
     *
     * === Case 2 === we are operating on a rotated lattice (Marked by 'C' and 'F')
     *
     *     current lattice              coarse lattice
     *   C--+--C--+--C--+--C          C-----C-----C-----C
     *   !  !  !  !  !  !  !          !     !     !     ! 
     *   +--F--+--F--+--F--+          !     !     !     ! 
     *   !  !  !  !  !  !  !          !     !     !     ! 
     *   C--+--C--+--C--+--C          C-----C-----C-----C
     *   !  !  !  !  !  !  !          !     !     !     !
     *   +--F--+--F--+--F--+          !     !     !     !
     *   !  !  !  !  !  !  !          !     !     !     !
     *   C--+--C--+--C--+--C          C-----C-----C-----C
     * 
     *  => Copy the C-points from the next coarser lattice (right)
     *
     * @param[in] phi_coarse Coarse sample to copy from
     * @param[in] phi_state Fine state to copy to (sample level as action)
     */
    virtual void copy_from_coarse(const std::shared_ptr<SampleState> phi_coarse,
                                  std::shared_ptr<SampleState> phi_state);

    /** @brief Copy coarse data points from state on finer level
     * 
     * As for the copy_from_coarse method, we need to consider two cases:
     * 
     * === Case 1 === we are operating on an un-rotated lattice
     * 
     *     current lattice               fine lattice
     *   C-----C-----C-----C          C--+--C--+--C--+--C          
     *   !     !     !     !          !  !  !  !  !  !  !          
     *   !     !     !     !          +--F--+--F--+--F--+          
     *   !     !     !     !          !  !  !  !  !  !  !          
     *   C-----C-----C-----C          C--+--C--+--C--+--C          
     *   !     !     !     !          !  !  !  !  !  !  !          
     *   !     !     !     !          +--F--+--F--+--F--+          
     *   !     !     !     !          !  !  !  !  !  !  !          
     *   C-----C-----C-----C          C--+--C--+--C--+--C          
     * 
     *  => Copy the C-points from the next finer lattice (right)
     * 
     * === Case 2 === we are operating on a rotated lattice
     *
     *     current lattice               fine lattice
     *   C--+--C--+--C--+--C          C--F--C--F--C--F--C
     *   !  !  !  !  !  !  !          !  !  !  !  !  !  !          
     *   +--C--+--C--+--C--+          F--C--F--C--F--C--F          
     *   !  !  !  !  !  !  !          !  !  !  !  !  !  !
     *   C--+--C--+--C--+--C          C--F--C--F--C--F--C          
     *   !  !  !  !  !  !  !          !  !  !  !  !  !  !
     *   +--C--+--C--+--C--+          F--C--F--C--F--C--F          
     *   !  !  !  !  !  !  !          !  !  !  !  !  !  !
     *   C--+--C--+--C--+--C          C--F--C--F--C--F--C          
     *      
     *  => Copy the C-points from the next finer lattice (right)
     * 
     * (note that the F-points will then later be filled in based on the
     * conditioned action)
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
    /** @brief Draw local value of state from heat bath at vertex (i,j)
     *
     * Update the local entry at position \f$\ell\f$ of the state using a heat bath defined by the neighbouring sites.
     * For this, observe that if all spins are kept fixed, except for the spin \f$\sigma_n\f$ at site \f$n\f$, then the
     * action can be written as:
     *
     * \f[
     * \begin{aligned}
     *   S_n &= -\beta \sigma_n\cdot \Delta_n\\
     *      &= -\beta |\Delta_n| \cos(\omega_n)\\
     *      &= \text{const.} + 2\beta |\Delta_n| \sin^2(\omega_n/2)
     * \end{aligned}
     * \f]
     *
     * Here \f$\omega_n\f$ is the angle between \f$\sigma_n\f$ and the sum of all neighbouring spins
     * \f$\Delta_n\f$. Hence, for the heat-bath update we need to proceed as follows
     *
     * 1. set \f$\sigma_n \mapsto \hat{\Delta}_n := |\Delta_n|^{-1}\Delta_n\f$
     * 2. draw an angle \f$\omega_n\f$ from the distribution \f$\pi_n\f$  with
     *   \f$\pi_n(\omega) \propto \exp\left[2\beta|\Delta_n| \sin^2(\omega/2)\right]\f$
     * 3. find a vector \f$\Delta_n^\perp\f$
     * 4. rotate \f$\sigma_n\f$ around \f$\Delta_n^\perp\f$ by the angle \f$\omega_n\f$
     * 5. draw an angle \f$\tau_n\f$ uniformly from the interval \f$[0,2\pi[\f$
     * 6. rotate \f$\sigma_n\f$ around \f$\Delta_n\f$ by an angle \f$\tau_n\f$
     *
     * Note that the final two steps leave the action invariant, and are hence also implemented in the
     * overrelaxation step.
     *
     *  @param[inout] phi_state State to update
     *  @param[in] i temporal index of site to update
     *  @param[in] j spatial index of site to update
     */
    void heatbath_ij_update(std::shared_ptr<SampleState> phi_state,
                            const int i,
                            const int j);


    /** @brief Extract 3d-vector at lattice vertex (i,j) from sample state and add it to vector
     *
     * @param[in] phi_state State \f$\phi\f$
     * @param[in] i temporal index
     * @param[in] j spatial index
     * @param[in] diag Use rotated lattice?
     * @param[out] sigma unit vector
     */
    void add_sigma(const std::shared_ptr<SampleState> phi_state,
                   const int i,
                   const int j,
                   Eigen::Vector3d& sigma) const {        
        double theta, phi;
        get_dofs(phi_state,i,j,theta,phi);
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
        if (rotated) {
            add_sigma(phi_state,i+1,j+1,Delta_n);
            add_sigma(phi_state,i+1,j-1,Delta_n);
            add_sigma(phi_state,i-1,j+1,Delta_n);
            add_sigma(phi_state,i-1,j-1,Delta_n);
        } else {
            add_sigma(phi_state,i+1,j  ,Delta_n);
            add_sigma(phi_state,i-1,j  ,Delta_n);
            add_sigma(phi_state,i  ,j+1,Delta_n);
            add_sigma(phi_state,i  ,j-1,Delta_n);
        }
        return Delta_n;
    }

protected:
    /** @brief Dimensionless coupling constant \f$\beta=1/(a_t a_x g^2)\f$*/
    const double beta;
    /** @brief Formulated on rotated lattice? */
    const bool rotated;
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
