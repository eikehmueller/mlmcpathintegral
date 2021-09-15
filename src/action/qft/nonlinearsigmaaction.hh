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
#include "action/clusteraction.hh"
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
 *      S[\sigma]= - \frac{1}{2} \beta \sum_n \sigma_n \cdot \Delta_n
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

class NonlinearSigmaAction : public QFTAction, public ClusterAction {
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
     */
    NonlinearSigmaAction(const std::shared_ptr<Lattice2D> lattice_,
                         const std::shared_ptr<Lattice2D> fine_lattice_,
                         const RenormalisationType renormalisation_,
                         const double beta_)
        : QFTAction(lattice_,fine_lattice_,renormalisation_),
          ClusterAction(lattice_),
          beta(beta_),
          uniform_dist(0.,2.*M_PI),
          neighbour_vertices(lattice_->get_neighbour_vertices()) {
              engine.seed(2481317);
              CoarseningType coarsening_type = lattice->get_coarsening_type();
              if (not (coarsening_type == CoarsenRotate) ) {
                  mpi_parallel::cerr << "ERROR: invalid coarsening for Non-linear sigma model." << std::endl;
                  mpi_parallel::cerr << "Has to be 'rotate'." << std::endl;
                  mpi_exit(EXIT_FAILURE);
                  throw std::runtime_error("...");
              }
          }

    /** @brief Return coupling constant \f$beta\f$ */
    double getbeta() const {
        return beta;
    }
            
    /** @brief return size of a sample, i.e. the number of vertices on the lattice */
    virtual unsigned int sample_size() const {
        return 2*lattice->getNvertices();
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
        std::shared_ptr<Lattice2D> coarse_lattice = lattice->get_coarse_lattice();       
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

    
    /** @brief Draw local value of state from heat bath at vertex with index ell
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
     *  @param[in] ell index of site to update
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

    /** @brief Extract and return 3d-vector at lattice vertex with index ell from state
     *
     * @param[in] phi_state State \f$\phi\f$
     * @param[in] ell linear index of vertex
     */
    Eigen::Vector3d get_sigma(const std::shared_ptr<SampleState> phi_state,
                              const unsigned int ell) const {        
        double theta=phi_state->data[2*ell];
        double phi=phi_state->data[2*ell+1];
        Eigen::Vector3d sigma = { sin(theta)*cos(phi),
                                  sin(theta)*sin(phi),
                                  cos(theta) };
        return sigma;
    }

    /** @brief Extract 3d-vector at lattice vertex with index ell from sample
     *         state and add it to vector. NB: the only difference to the set_sigma()
     *         is that here the vector is inremented, not just set.
     *
     * @param[in] phi_state State \f$\phi\f$
     * @param[in] ell linear index of vertex
     * @param[out] sigma unit vector
     */
    void add_sigma(const std::shared_ptr<SampleState> phi_state,
                   const unsigned int ell,
                   Eigen::Vector3d& sigma) const {        
        double theta=phi_state->data[2*ell];
        double phi=phi_state->data[2*ell+1];
        sigma[0] += sin(theta)*cos(phi);
        sigma[1] += sin(theta)*sin(phi);
        sigma[2] += cos(theta);
    }

    /** @brief Extract 3d-vector which is given by the sum of the 3d field vectors of
     * the neighbours of site with linex vertex index ell
     *
     * @param[in] phi_state State \f$\phi\f$
     * @param[in] ell linear index of vertex
     * @param[out] Delta_n resulting vector
     */
    Eigen::Vector3d delta_neighbours(const std::shared_ptr<SampleState> phi_state,
                                     const unsigned int ell) const {
        Eigen::Vector3d Delta_n(0.,0.,0.);
        for (int k=0;k<4;++k) {
            add_sigma(phi_state,neighbour_vertices[ell][k],Delta_n);
        }
        return Delta_n;
    }
        
    /** @brief Change \f$S_{\ell}\f$ in energy used in bonding probabilities
     *
     * The probability to have a bond between sites \f$i\f$ and \f$j\f$
     * is given by \f$1-e^{\min(0,-S_{\ell})}\f$
     *
     * @param[in] x_path Sample state
     * @param[in] i index of first vertex
     * @param[in] i index of second vertex
     */
    virtual double S_ell(const std::shared_ptr<SampleState> phi_state,
                         const unsigned int i,
                         const unsigned int j) const;

    
    /** @brief Flip the spin at a given site with respect to the reflection plane
     *
     * @param[inout] phi_state State to process
     * @param[in] ell Vertex at which to flip the spin
     */
    virtual void flip(std::shared_ptr<SampleState> phi_state,
                      const unsigned int ell) const;

    /** @brief Draw new reflection angle (normal for reflection plan) for cluster sampler */
    virtual void new_reflection() const;

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
    /** @brief Reference to the neighbour-list of the underlying lattice */
    const std::vector<std::vector<unsigned int> >& neighbour_vertices;
    /** @brief spin-flip vector for cluster updates */
    mutable Eigen::Vector3d sigma_spinflip;
};

#endif // NONLINEARSIGMAACTION_HH
