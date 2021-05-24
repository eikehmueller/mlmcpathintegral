#ifndef QUENCHEDSCHWINGERCONDITIONEDFINEACTION_HH
#define QUENCHEDSCHWINGERCONDITIONEDFINEACTION_HH
#include <random>
#include <memory>
#include <algorithm>
#include "common/auxilliary.hh"
#include "common/samplestate.hh"
#include "mpi/mpi_random.hh"
#include "lattice/lattice2d.hh"
#include "action/action.hh"
#include "action/conditionedfineaction.hh"
#include "action/qft/quenchedschwingeraction.hh"
#include "distribution/expcosdistribution.hh"
#include "distribution/besselproductdistribution.hh"
#include "distribution/approximatebesselproductdistribution.hh"
#include "distribution/gaussianfillindistribution.hh"

/** @file quenchedschwingerconditionedfineaction.hh
 * @brief Header file for quenched Schwinger model conditioned fine action class
 */

/** @class QuenchedSchwingerConditionedFineAction
 *
 * @brief Conditioned fine action for the quenched Schwinger model
 */
class QuenchedSchwingerConditionedFineAction : public ConditionedFineAction {
public:
    /** @brief Constructor
     *
     * Construct new instance
     *
     * @param[in] action_ Underlying action class
     */
    QuenchedSchwingerConditionedFineAction (const std::shared_ptr<QuenchedSchwingerAction> action_) :
        action(action_),
        beta(action->getbeta()),
        uniform_dist(-M_PI,+M_PI),
        exp_cos_dist(action->getbeta()),
        bessel_product_dist(NULL),
        gaussian_fillin_dist(NULL),
        approximate_bessel_product_dist(NULL) {
            if (beta>8.0) {
                // Use Gaussian approximation for larger values of beta?
                bool gaussian_approximation = false;
                if (gaussian_approximation) {
                    bool add_gaussian_noise = true;
                    gaussian_fillin_dist = std::make_shared<GaussianFillinDistribution>(beta,
                                                                                        add_gaussian_noise);                                                                                    
                } else {
                    approximate_bessel_product_dist = std::make_shared<ApproximateBesselProductDistribution>(beta);
                }
            } else {
                bessel_product_dist = std::make_shared<BesselProductDistribution>(beta);
            }
            // Sanity check: make sure at least one fill-in distribution is defined
            if ( (gaussian_fillin_dist == NULL) and 
                 (approximate_bessel_product_dist == NULL) and 
                 (bessel_product_dist == NULL) ) {
                mpi_parallel::cerr << "ERROR: no fill-in distribution defined!" << std::endl;
                mpi_exit(EXIT_FAILURE);
            }
            engine.seed(71814151);
        }

    /** @brief Destructor */
    virtual ~QuenchedSchwingerConditionedFineAction() {}

    /** @brief Fill in fine points
     *
     * Given a path \f$\phi\f$ for which the coarse links have been set, fill in
     * the fine links by sampling from the conditioned action.
     *
     * @param[in] phi_state_n State \f$\phi_n\f$ at previous step (for pCN)
     * @param[inout] phi_state State \f$\phi\f$ to fill
     */
    virtual void fill_fine_points(const std::shared_ptr<SampleState> phi_state_n,
                                  std::shared_ptr<SampleState> phi_state) const;

    /** @brief Evaluate conditioned action at fine points
     *
     * Given a path \f$\phi\f$ for which all points have been set, evaluate the
     * conditioned action.
     *
     * @param[inout] phi_state Path \f$\phi\f$ to evaluate
     */
    virtual double evaluate(const std::shared_ptr<SampleState> phi_state) const;

private:
    /** @brief Underlying action class */
    const std::shared_ptr<QuenchedSchwingerAction> action;
    /** @brief Coupling constant \f$\beta\f$ */
    const double beta;
    /** @brief Random number engine */
    typedef mpi_parallel::mt19937_64 Engine;
    /** @brief Type of Mersenne twister engine */
    mutable Engine engine;
    /** @brief Uniform probability distribution for drawing links on perimeter */
    mutable std::uniform_real_distribution<double> uniform_dist;
    /** @brief Probability distribution for drawing horizontal interior links */
    mutable ExpCosDistribution exp_cos_dist;
    /** @brief Probability distribution for drawing vertical interior links */
    mutable std::shared_ptr<BesselProductDistribution> bessel_product_dist;
    /** @brief Probability distribution for drawing all using pCN */
    mutable std::shared_ptr<GaussianFillinDistribution> gaussian_fillin_dist;
    /** @brief Approximate probability distribution for drawing vertical interior links */
    mutable std::shared_ptr<ApproximateBesselProductDistribution> approximate_bessel_product_dist;

};

struct QuenchedSchwingerConditionedFineActionFactory : public ConditionedFineActionFactory {
    /** @brief Destructor */
    virtual ~QuenchedSchwingerConditionedFineActionFactory() {}

    /** @brief Method for constructing a conditioned fine action
     * @param[in] action Coarse level action
     */
    virtual std::shared_ptr<ConditionedFineAction> get(std::shared_ptr<Action> action) {
        return std::make_shared<QuenchedSchwingerConditionedFineAction>(std::dynamic_pointer_cast<QuenchedSchwingerAction>(action));
    }
};


#endif // QUENCHEDSCHWINGERCONDITIONEDFINEACTION_HH
