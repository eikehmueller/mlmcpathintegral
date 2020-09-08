#ifndef SAMPLER_HH
#define SAMPLER_HH SAMPLER_HH
#include <memory>
#include <iostream>
#include <iomanip>
#include <vector>
#include "common/samplestate.hh"
#include "montecarlo/mcmcstep.hh"
#include "action/action.hh"

/** @file sampler.hh
 * @brief Header file for sampler base class
 */

/** @class Sampler
 * @brief Base class for sampler
 *
 * Abstract base class for sampler from distribution over states
 */
class Sampler : public MCMCStep {
public:
  /** @brief Create new instance
   *
   */
  Sampler() : MCMCStep() {}

  /** @brief Delete instance */
  virtual ~Sampler() {};

  /** @brief Draw a sample 
   *
   * returns a sample state \f$\phi\f$
   *
   * @param[out] phi_state State \f$\phi\f$ drawn from distribution
   */
  virtual void draw(std::shared_ptr<SampleState> phi_state) = 0;
};

class SamplerFactory {
public:
  /** @brief extract a sampler for a given action */
  virtual std::shared_ptr<Sampler> get(std::shared_ptr<Action> action) = 0;
};


#endif // SAMPLER_HH
