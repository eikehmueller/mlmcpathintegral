#ifndef FASTBESSEL_HH
#define FASTBESSEL_HH FASTBESSEL_HH

#include <cmath>
#include <gsl/gsl_sf_bessel.h>

/** @file fastbessel.hh
 *
 * @brief Fast implementation of modified Bessel function \f$I_0\f$
 * for large arguments
 */

/** @brief Evaluate Bessel function \f$I_0(z)\f$ scaled by \f$e^{-z}\f$
 *
 * This uses an appropriate GSL function for small arguments and an
 * asymptotic expansion for large arguments \f$z>z_0\f$. More specifically,
 * the following (truncation) series representation is used:
 *
 * \f[
 *   e^{-z} I_0(z) \approx \frac{1}{2\pi z_{\text{inv}}} \sum_{k=0}^{k_{\max}}
 * a_k z^{-k} \f]
 *
 * The coefficients are given by \f$a_k=\frac{((2k-1)!!)^2}{8^k k!}\f$.
 */

/** @brief Evaluate scaled modified Bessel function \f$e^{-z}I_0(z)\f$
 *
 * @param[in] z Argument at which to evaluate \f$e^{-z}I_0\f$
 */
const double fast_bessel_I0_scaled(const double z);

/** @brief Template for evaluating the asymptotic expansion coefficients
 *
 * The asymptotic expansion coefficient \f$a_k\f$ satisfy the recursion relation
 * \f$a_k = \frac{(2k-1)^2}{8n}a_{n-1}\f$ with \f$a_0\f$ This is implemented
 * using template metaprogramming.
 */
template <int N> struct ModifiedBesselCoefficient {
  static constexpr double value = 0.125 * (2.0 * N - 1.0) * (2.0 * N - 1.0) /
                                  N * ModifiedBesselCoefficient<N - 1>::value;
};

/** @brief Base case template for evaluating the asymptotic expansion
 * coefficient \f$a_0\f$
 *
 * The base case is given by \f$a_0=1\f$.
 */
template <> struct ModifiedBesselCoefficient<0> {
  static constexpr double value = 1.0;
};

#endif // FASTBESSEL_HH
