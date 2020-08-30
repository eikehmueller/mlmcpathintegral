#ifndef SAMPLER_HH
#define SAMPLER_HH SAMPLER_HH
#include <memory>
#include <iostream>
#include <iomanip>
#include <vector>
#include "montecarlo/mcmcstep.hh"
#include "action/action.hh"
#include "fields/path.hh"

/** @file sampler.hh
 * @brief Header file for sampler base class
 */

/** @class Sampler
 * @brief Base class for sampler
 *
 * Abstract base class for sampler from distribution over paths
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
   * returns a sample path \f$X\f$
   *
   * @param[out] x_path Path \f$X\f$ drawn from distribution
   */
  virtual void draw(std::shared_ptr<Path> x_path) = 0;
};

class SamplerFactory {
public:
  /** @brief extract a sampler for a given action */
  virtual std::shared_ptr<Sampler> get(std::shared_ptr<Action> action) = 0;
};


#endif // SAMPLER_HH
