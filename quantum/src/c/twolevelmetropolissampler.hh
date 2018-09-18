#ifndef TWOLEVELMETROPOLISSAMPLER_HH
#define TWOLEVELMETROPOLISSAMPLER_HH TWOLEVELMETROPOLISSAMPLER_HH
#include <memory>
#include <random>
#include "action.hh"
#include "mcmcstep.hh"
#include "conditionedfineaction.hh"

/** @file twolevelmetropolissampler.hh
 * @brief Header file for two level Metropolis sampler
 */

/** @class TwoLevelMetropolisSampler 
 *
 * @brief Two level sampler
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
class TwoLevelMetropolisSampler : public MCMCStep {
public:
  /** @brief Create a new instance
   *
   * @param[in] coarse_sampler_ Sampler on coarse level \f$\ell-1\f$
   * @param[in] coarse_action_ Action on coarse level \f$\ell-1\f$
   * @param[in] fine_action_ Action on fine level \f$\ell\f$
   * @param[in] conditioned_fine_action_ Conditioned fine action object for
   *            filling in the fine points
   */
  TwoLevelMetropolisSampler(std::shared_ptr<Sampler> coarse_sampler_,
                            const std::shared_ptr<Action> coarse_action_,
                            const std::shared_ptr<Action> fine_action_,
                            const std::shared_ptr<ConditionedFineAction> conditioned_fine_action_) :
    MCMCStep(),
    coarse_sampler(coarse_sampler_),
    coarse_action(coarse_action_),
    fine_action(fine_action_),
    conditioned_fine_action(conditioned_fine_action_) {
    assert(2*coarse_action->getM_lat()==fine_action->getM_lat());
    theta_coarse = std::make_shared<Path>(coarse_action->getM_lat(),
                                          coarse_action->getT_final());
    theta_fine = std::make_shared<Path>(fine_action->getM_lat(),
                                        coarse_action->getT_final());
    theta_fine_C = std::make_shared<Path>(coarse_action->getM_lat(),
                                          coarse_action->getT_final());
    theta_prime = std::make_shared<Path>(fine_action->getM_lat(),
                                         coarse_action->getT_final());
    engine.seed(89216491);
    reset_stats();
  }

  /** @brief draw new fine-/coarse-level path pair
   *
   * @param[out] x_path Vector of pointers to fine- and coarse- path
   */
  virtual void draw(std::vector<std::shared_ptr<Path>> x_path);
                              
protected:
  /** @brief Sampler on coarse level */
  std::shared_ptr<Sampler> coarse_sampler;
  /** @brief Action on coarse level */
  const std::shared_ptr<Action> coarse_action;
  /** @brief Action on fine level */
  const std::shared_ptr<Action> fine_action;
  /** @brief Conditioned fine action */
  const std::shared_ptr<ConditionedFineAction> conditioned_fine_action;
  /** @brief Temporary state vector on coarse level \f$\theta^n_{\ell-1}\f$ */
  mutable std::shared_ptr<Path> theta_coarse;
  /** @brief Temporary state vector on fine level \f$\theta^n_{\ell}\f$ */
  mutable std::shared_ptr<Path> theta_fine;
  /** @brief Coarse part of state vector on fine level \f$\theta^n_{\ell,C}\f$ */
  mutable std::shared_ptr<Path> theta_fine_C;
  /** @brief Trial state vector on fine level \f$\theta'_{\ell}\f$ */
  mutable std::shared_ptr<Path> theta_prime;
  /** @brief Random number engine */
  typedef std::mt19937_64 Engine;
  /** @brief Type of Mersenne twister engine */
  mutable Engine engine;
  /** @brief Type of uniform distribution */
  typedef std::uniform_real_distribution<double> Uniform;
  /** @brief Uniform distribution in [0,1] for accept/reject step */
  mutable Uniform uniform_dist;
};
#endif // TWOLEVELMETROPOLISSAMPLER_HH
