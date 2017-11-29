#ifndef SAMPLER_HH
#define SAMPLER_HH SAMPLER_HH
#include "path.hh"
#include <vector>

/** @file sampler.hh
 * @brief Header file for sampler base class
 */

/** @class Sampler
 * @brief Base class for sampler
 *
 * Abstract base class for sampler from distribution over paths
 */
class Sampler {
public:
  /** @brief Create new instance
   *
   * @param[in] M_lat_ Number \f$M\f$ of time slices
   */
  Sampler(const unsigned int M_lat_) : M_lat(M_lat_) {}

  /** @brief Return number of timeslices */
  unsigned int getM_lat() const { return M_lat; }

  /** @brief Draw a sample 
   *
   * returns a sample path \f$X\f$
   *
   * @param[out] x_path Path \f$X\f$ drawn from distribution
   */
  const virtual void draw(std::vector<Path*> x_path) = 0;

protected:
  /** @brief Number of time slices */
  const unsigned int M_lat;
};

#endif // SAMPLER_HH
