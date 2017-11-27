#ifndef SAMPLER_HH
#define SAMPLER_HH SAMPLER_HH
#include <vector>

/** @class Sampler
 * @brief Base class for sampler
 *
 * Abstract base class for sampler from distribution over paths
 */
class Sampler {
public:
  /** @brief Create new instance
   *
   * @param M_lat Number \f$M\f$ of time slices
   */
  Sampler(const unsigned int M_lat_) : M_lat(M_lat_) {}

  /** @brief Return number of timeslices */
  unsigned int getM_lat() const { return M_lat; }

  /** @brief Draw a sample 
   *
   * returns a sample path \f$X\f$
   * @param[out] x_path Path \f$X\f$ drawn from distribution
   */
  const virtual void draw(std::vector<Path*> X_path) = 0;

protected:
  const unsigned int M_lat;
};

#endif // SAMPLER_HH
