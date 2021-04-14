#ifndef GAUSSIANCONDITIONEDFINEACTION_HH
#define GAUSSIANCONDITIONEDFINEACTION_HH GAUSSIANCONDITIONEDFINEACTION_HH
#include <random>
#include <memory>
#include <algorithm>
#include "common/auxilliary.hh"
#include "mpi/mpi_random.hh"
#include "common/samplestate.hh"
#include "action/qm/qmaction.hh"
#include "action/conditionedfineaction.hh"

/** @file gaussianconditionedfineaction.hh
 * @brief Header file for Gaussian conditioned fine action class
 */

/** @class GaussianConditionedFineAction
 *
 * @brief Gaussian conditioned fine action
 *
 * The conditioned action is given by
 * \f[
 *    S^{cond}(X;X_-,X_+) = \frac{1}{2}W''(X_-,X_+)\left(X-\overline{X}(X_-,X_+)\right)^2
 * \f]
 * where the curvature \f$W''(X_-,X_+)\f$ and mean \f$\overline{X}(X_-,X_+)\f$
 * are given by the getWcurvature() and getWminimum() methods of underlying
 * action class which is passed to the constructor.
 */
class GaussianConditionedFineAction : public ConditionedFineAction {
public:
    /** @brief Constructor
     *
     * Construct new instance
     *
     * @param[in] action_ Underlying action class
     */
    GaussianConditionedFineAction (const std::shared_ptr<QMAction> action_) : action(action_) {
        engine.seed(11897197);
    }

    /** @brief Destructor */
    virtual ~GaussianConditionedFineAction() {}

    /** @brief Fill in fine points
     *
     * Given a path \f$X\f$ for which the coarse points have been set, fill in
     * the fine points by sampling from the conditioned action.
     *
     * @param[in] x_path_n SampleState \f$X_n\f$ at previous step
     * @param[inout] x_path SampleState \f$X\f$ to fill
     */
    virtual void fill_fine_points(const std::shared_ptr<SampleState> x_path_n,
                                  std::shared_ptr<SampleState> x_path) const;

    /** @brief Evaluate conditioned action at fine points
     *
     * Given a path \f$X\f$ for which all points have been set, evaluate the
     * conditioned action.
     *
     * @param[inout] x_path SampleState \f$X\f$ to fill
     */
    virtual double evaluate(const std::shared_ptr<SampleState> x_path) const;

private:
    /** @brief Underlying action class */
    const std::shared_ptr<QMAction> action;
    /** @brief Random number engine */
    typedef mpi_parallel::mt19937_64 Engine;
    /** @brief Type of Mersenne twister engine */
    mutable Engine engine;
    /** @brief Type of normal distribution */
    typedef std::normal_distribution<double> Normal;
    /** @brief Normal distribution for drawing from distribution */
    mutable Normal normal_dist;

};

struct GaussianConditionedFineActionFactory : public ConditionedFineActionFactory {
    /** @brief Destructor */
    virtual ~GaussianConditionedFineActionFactory() {}

    /** @brief Method for constructing a conditioned fine action
     * @param[in] action Coarse level action
     */
    virtual std::shared_ptr<ConditionedFineAction> get(std::shared_ptr<Action> action) {
        return std::make_shared<GaussianConditionedFineAction>(std::dynamic_pointer_cast<QMAction>(action));
    }
};


#endif // CONDITIONEDFINEACTION_HH
