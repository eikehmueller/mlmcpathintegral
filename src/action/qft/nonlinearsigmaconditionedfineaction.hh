#ifndef NONLINEARSIGMACONDITIONEDFINEACTION_HH
#define NONLINEARSIGMACONDITIONEDFINEACTION_HH NONLINEARSIGMACONDITIONEDFINEACTION_HH
#include <random>
#include <memory>
#include <algorithm>
#include "common/auxilliary.hh"
#include "common/samplestate.hh"
#include "lattice/lattice2d.hh"
#include "distribution/expsin2distribution.hh"
#include "action/action.hh"
#include "action/conditionedfineaction.hh"
#include "action/qft/nonlinearsigmaaction.hh"

/** @file nonlinearsigmaconditionedfineaction.hh
 * @brief Header file for non-linear sigma model conditioned fine action class
 */

/** @class NonlinearSigmaConditionedFineAction
 *
 * @brief Conditioned fine action for the non-linear sigma model
 * 
 * Fills in the values at the sites marked by 'X', given the values at the coarse
 * level sites marked by 'C' in the following diagrams
 * 
 *    unrotated action              rotated action
 *   C--X--C--X--C--X--C         C--+--C--+--C--+--C
 *   !  !  !  !  !  !  !         !  !  !  !  !  !  !
 *   X--C--X--C--X--C--X         +--X--+--X--+--X--+
 *   !  !  !  !  !  !  !         !  !  !  !  !  !  !
 *   C--X--C--X--C--X--C         C--+--C--+--C--+--C
 *   !  !  !  !  !  !  !         !  !  !  !  !  !  !
 *   X--C--X--C--X--C--X         +--X--+--X--+--X--+
 *   !  !  !  !  !  !  !         !  !  !  !  !  !  !
 *   C--X--C--X--C--X--C         C--+--C--+--C--+--C
 * 
 * Note that in the unrotated case the sites that need to be filled in
 * have exactly one even and one odd index, i.e. i+j is odd. In the rotated
 * case, but i *and* j are odd.
 * 
 */
class NonlinearSigmaConditionedFineAction : public ConditionedFineAction {
public:
    /** @brief Constructor
     *
     * Construct new instance
     *
     * @param[in] action_ Underlying action class
     */
    NonlinearSigmaConditionedFineAction (const std::shared_ptr<NonlinearSigmaAction> action_) :
        action(action_),
        beta(action->getbeta()) { engine.seed(115147); }

    /** @brief Destructor */
    virtual ~NonlinearSigmaConditionedFineAction() {}

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
    const std::shared_ptr<NonlinearSigmaAction> action;
    /** @brief Coupling constant \f$\beta\f$ */
    const double beta;
    /** @brief Random number engine */
    mutable mpi_parallel::mt19937_64 engine;
    /** @brief Probability distribution for evaluating conditional probabilities */
    mutable ExpSin2Distribution exp_sin2_dist;
};

struct NonlinearSigmaConditionedFineActionFactory : public ConditionedFineActionFactory {
    /** @brief Destructor */
    virtual ~NonlinearSigmaConditionedFineActionFactory() {}

    /** @brief Method for constructing a conditioned fine action
     * @param[in] action Coarse level action
     */
    virtual std::shared_ptr<ConditionedFineAction> get(std::shared_ptr<Action> action) {
        std::shared_ptr<NonlinearSigmaAction> nonlinearsigma_action = std::dynamic_pointer_cast<NonlinearSigmaAction>(action);
        // Select correct type based on coarsening type
        return std::make_shared<NonlinearSigmaConditionedFineAction>(nonlinearsigma_action);
    }
};

#endif // NONLINEARSIGMACONDITIONEDFINEACTION_HH
