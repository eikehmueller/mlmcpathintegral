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
        bessel_product_dist(action->getbeta()) {
        engine.seed(71814151);
    }

    /** @brief Destructor */
    virtual ~QuenchedSchwingerConditionedFineAction() {}

    /** @brief Fill in fine points
     *
     * Given a path \f$\phi\f$ for which the coarse links have been set, fill in
     * the fine links by sampling from the conditioned action.
     *
     * @param[inout] phi_state State \f$\phi\f$ to fill
     */
    virtual void fill_fine_points(std::shared_ptr<SampleState> phi_state) const;

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
    mutable BesselProductDistribution bessel_product_dist;
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
