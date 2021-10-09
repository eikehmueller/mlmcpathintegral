#ifndef GFFCONDITIONEDFINEACTION_HH
#define GFFCONDITIONEDFINEACTION_HH GFFCONDITIONEDFINEACTION_HH
#include <random>
#include <memory>
#include <algorithm>
#include "common/samplestate.hh"
#include "lattice/lattice2d.hh"
#include "action/action.hh"
#include "action/qft/gffaction.hh"
#include "action/conditionedfineaction.hh"

/** @file gffconditionedfineaction.hh
 * @brief Header file for 2d Gaussian Free Field (GFF)
 */

/** @class GFFConditionedFineAction
 *
 * @brief Conditioned fine action for the Gaussian Free Field (GFF)
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
class GFFConditionedFineAction : public ConditionedFineAction {
public:
    /** @brief Constructor
     *
     * Construct new instance
     *
     * @param[in] action_ Underlying action class
     */
    GFFConditionedFineAction (const std::shared_ptr<GFFAction> action_) :
        action(action_),
        mu2(action_->mu2) { engine.seed(115147); }

    /** @brief Destructor */
    virtual ~GFFConditionedFineAction() {}

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
    const std::shared_ptr<GFFAction> action;
    /** @brief Random number engine */
    mutable mpi_parallel::mt19937_64 engine;
    /** @brief Normal distribution used for sampling */
    mutable std::normal_distribution<double> normal_dist;
    /** @brief Squared mass parameter \f$\mu^2\f$ */
    const double mu2;
};

struct GFFConditionedFineActionFactory : public ConditionedFineActionFactory {
    /** @brief Destructor */
    virtual ~GFFConditionedFineActionFactory() {}

    /** @brief Method for constructing a conditioned fine action
     * @param[in] action Coarse level action
     */
    virtual std::shared_ptr<ConditionedFineAction> get(std::shared_ptr<Action> action) {
        std::shared_ptr<GFFAction> gff_action = std::dynamic_pointer_cast<GFFAction>(action);
        // Select correct type based on coarsening type
        return std::make_shared<GFFConditionedFineAction>(gff_action);
    }
};

#endif // GFFCONDITIONEDFINEACTION_HH
