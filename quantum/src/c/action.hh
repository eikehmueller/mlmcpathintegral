#ifndef ACTION_HH
#define ACTION_HH ACTION_HH
#include <cassert>
#include <random>
#include <vector>
#include <Eigen/Dense>
#include "path.hh"
#include "sampler.hh"

/** @class Action
 *
 * @brief Base class for action
 *
 * Allows calculation of action
 *   \f$S[X]=\sum_{j=0}^{M-1}\left(\frac{m_0}{2}\frac{(X_{j+1}-X_j)}{a^2}+V(x)\right)\f$ 
 * for one-dimensional quantum problem with periodic boundary conditions.
 */
class Action {
public:
  /* @brief Initialise class
   *
   * 
   * @param[in] M_lat Number of time slices \f$M\f$
   * @param[in] T_final Final time \f$T\f$
   * @param[in] m0: Mass of particle \f$m_0\f$
   */
  Action(const unsigned int M_lat_,
         const double T_final_,
         const double m0_)
    : M_lat(M_lat_), T_final(T_final_), m0(m0_), a_lat(T_final_/M_lat_) {
    assert(m0>0.0);
    assert(T_final>0.0);
  }

  /** @brief Return number of timeslices */
  unsigned int getM_lat() const { return M_lat; }

  /** @brief Return mass \f$m_0\f$ */
  double getm0() const { return m0;}

  /** @brief Return lattice spacing \f$a\f$ */
  double geta_lat() const { return a_lat;}

  /* @brief Evaluate action for a specific path
   * 
   * Calculate \$S[X]\f$ for a specific path
   *
   * @param x_path: Path, has to be am array of length \f$M\f$
   */
  const double virtual evaluate(const Path* x_path) const = 0;
protected:
  const unsigned int M_lat;
  const double T_final;
  const double m0;
  const double a_lat;
};

#endif // ACTION_HH
