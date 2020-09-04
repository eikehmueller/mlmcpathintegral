#ifndef QOISUSCEPTIBILITY_HH
#define QOISUSCEPTIBILITY_HH QOISUSCEPTIBILITY_HH
#include <memory>
#include <cmath>
#include "common/auxilliary.hh"
#include "common/samplestate.hh"
#include "lattice/lattice1d.hh"
#include "qoi/quantityofinterest.hh"

/** @file qoisusceptibility.hh
 * @brief Header file for QoI Q^2/T where Q is the topological charge
 */

/** @class QoISusceptibility
 *
 * @brief class for calculating the susceptibility for the QM rotor
 *
 * For a path \f$X\f$ the susceptibility is defined as
 * \f$\chi_t = Q[X]^2/T \f$ where
 * \f[
 *   Q[X] = \frac{1}{2\pi} \sum_{j=0}^{M_{lat}-1} (X_j - X_{j-1}) \mod [-\pi,\pi]
 * \f]
 * and the time is \f$T=a_{lat}M_{lat}\f$.
 */

class QoISusceptibility : public QoI {
public: 
  /** @brief Create new instance 
   *
   * @param[in] M_lat_ Number of time slices \f$M\f$
   */
  QoISusceptibility(const std::shared_ptr<Lattice1D> lattice) : 
    four_pi2_inv(0.25/(M_PI*M_PI)),
    T_final(lattice->getT_final()), 
    M_lat(lattice->getM_lat()) {}

  /** @brief Destructor */
  virtual ~QoISusceptibility() {}

  /** @brief Evaluate on a path
   *
   * @param[in] x_path Path \f$X\f$ on which to evaluate the QoI
   */
  const double virtual evaluate(const std::shared_ptr<SampleState> x_path);
  
private:
  /** @brief Constant \f$1/(4\pi^2)\f$ */
  const double four_pi2_inv;
  /** @brief Physical lattice size */
  const double T_final;
  /** @brief Number of lattice points */
  const unsigned int M_lat;
};

#endif // QOISUSCEPTIBILITY_HH
