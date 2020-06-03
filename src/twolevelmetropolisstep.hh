#ifndef TWOLEVELMETROPOLISSTEP_HH
#define TWOLEVELMETROPOLISSTEP_HH TWOLEVELMETROPOLISSTEP_HH
#include <memory>
#include <random>
#include "action.hh"
#include "mcmcstep.hh"
#include "conditionedfineaction.hh"
#include "mpi_random.hh"
#include "timer.hh"

/** @file twolevelmetropolisstep.hh
 * @brief Header file for two level Metropolis step
 */

/** @class TwoLevelMetropolisStep
 *
 * @brief Two level Metropolis step
 *
 * Given a way of creating coarse samples \f$\theta_{\ell-1}\f$ from a 
 * sampler on level \f$\ell-1\f$, create fine samples \f$\theta_{\ell}\f$
 * using a two level Metropolis Hastings process.
 * For this, states on level \f$\ell\f$ is written as 
 * \f$\theta_\ell = [\theta_{\ell,C},\theta_{\ell,F}]\f$ where
 * \f$\theta_{\ell,C}\f$ are the states on the even timeslices and 
 * \f$\theta_{\ell,F}\f$ are the states on the odd timeslices.
 *
 * To construct a new trial state, take the even modes from the next level
 * \f$\ell\f$ sample, \f$\theta'_{\ell,C}=\theta_{\ell}^{n+1}\f$. The 
 * odd modes are drawn from the free action, conditioned on the two
 * neighbouring odd modes \f$x_-\f$ and \f$x_+\f$, i.e. 
 *
 * \f[
   p(x) \sim \exp\left[-\frac{m_0}{a} \left(x-\frac{x_-+x_+}{2}\right)^2\right].
   \f]
 *
 * Finally, accept and reject according to \f$\min\left\{1,\beta\right\}\f$ with
 *
 * \f[
      \beta = \frac{\pi^{\ell}(\theta'_\ell)\pi^{\ell-1}(\theta^n_{\ell,C})\pi_{free}^{\ell}(\theta^n_{\ell,F}|\theta^n_{\ell,C})}{\pi^{\ell}(\theta^n_\ell)\pi^{\ell-1}(\theta'_{\ell,C})\pi_{free}^{\ell}(\theta'_{\ell,F}|\theta'_{\ell,C})}
 * \f]
 *
 * This guarantees that the fine level samples have the correct distribution.
 */
class TwoLevelMetropolisStep : public MCMCStep {
public:
  /** @brief Create a new instance
   *
   * @param[in] coarse_action_ Action on coarse level \f$\ell-1\f$
   * @param[in] fine_action_ Action on fine level \f$\ell\f$
   * @param[in] conditioned_fine_action_ Conditioned fine action object for
   *            filling in the fine points
   */
  TwoLevelMetropolisStep(const std::shared_ptr<Action> coarse_action_,
                         const std::shared_ptr<Action> fine_action_,
                         const std::shared_ptr<ConditionedFineAction> conditioned_fine_action_);

  /** @brief Destructor */
  virtual ~TwoLevelMetropolisStep() {}

  /** @brief draw new fine path given a coarse path
   *
   * @param[out] x_coarse_path Coarse path
   * @param[out] x_path Resulting fine path
   */
  virtual void draw(const std::shared_ptr<Path> x_coarse_path,
                    std::shared_ptr<Path> x_path);
               
  /** @brief Set current state to particular value
   *
   * @param[in] x_path
   */
  virtual void set_state(std::shared_ptr<Path> x_path);

  /** Return cost per sample */
  virtual double cost_per_sample() { return cost_per_sample_; }
  
protected:
  /** @brief Action on coarse level */
  const std::shared_ptr<Action> coarse_action;
  /** @brief Action on fine level */
  const std::shared_ptr<Action> fine_action;
  /** @brief Conditioned fine action */
  const std::shared_ptr<ConditionedFineAction> conditioned_fine_action;
  /** @brief Temporary state vector on fine level \f$\theta^n_{\ell}\f$ */
  mutable std::shared_ptr<Path> theta_fine;
  /** @brief Coarse part of state vector on fine level \f$\theta^n_{\ell,C}\f$ */
  mutable std::shared_ptr<Path> theta_fine_C;
  /** @brief Trial state vector on fine level \f$\theta'_{\ell}\f$ */
  mutable std::shared_ptr<Path> theta_prime;
  /** @brief Random number engine */
  typedef mpi_parallel::mt19937_64 Engine;
  /** @brief Type of Mersenne twister engine */
  mutable Engine engine;
  /** @brief Type of uniform distribution */
  typedef std::uniform_real_distribution<double> Uniform;
  /** @brief Uniform distribution in [0,1] for accept/reject step */
  mutable Uniform uniform_dist;
  /** @brief cost per sample */
  mutable double cost_per_sample_;
  /** @brief cached value of fine action evaluated for current state */
  mutable double fine_action_theta_fine;
  /** @brief cached value of conditioned fine action evaluated for current state */
  mutable double conditioned_fine_action_theta_fine;
};
#endif // TWOLEVELMETROPOLISSTEP_HH
