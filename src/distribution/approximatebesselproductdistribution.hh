#ifndef APPROXIMATEBESSELPRODUCTDISTRIBUTION_HH
#define APPROXIMATEBESSELPRODUCTDISTRIBUTION_HH                                \
  APPROXIMATEBESSELPRODUCTDISTRIBUTIONHH
#include "common/auxilliary.hh"
#include "common/fastbessel.hh"
#include "mpi/mpi_wrapper.hh"
#include <algorithm>
#include <cmath>
#include <iostream>
#include <random>

/** @file approximatebesselproductdistribution.hh
 * @brief Header file for approximation of BesselProductDistribution
 */

/** @class ApproximateBesselProductDistribution
 *
 * @brief Class for sampling from the distribution which approximates
 * BesselProductDistribution for large values of the coupling constant
 * \f$\beta\f$.
 *
 * To approximately sample from \f$p(x|x_+,x_-)\f$, this class allows sampling
 * from the distribution \f$\tilde{p}(x|x_+,x_-)\approx p(x|x_+,x_-)\f$, which
 * can be written as a sum of two Gaussians:
 *
 * \f[
 *  \tilde{p}(x|x_0,0) = N_+(x_0) S_+(x|x_0) + N_-(x_0) S_-(x|x_0)
 * \f]
 *
 * with
 *
 * \f[
 *    S_+(x|x_0) = \frac{1}{2\pi\sigma_+^2}\sum_{k\in\mathbb{Z}}
 * \exp\left[-\frac{1}{2\sigma_+^2}\left(x-\frac{x_0}{2}+2k\pi\right)^2\right]\\
 *    S_-(x|x_0) = \frac{1}{2\pi\sigma_-^2}\sum_{k\in\mathbb{Z}}
 * \exp\left[-\frac{1}{2\sigma_-^2}\left(x-\frac{x_0}{2}+\pi+2k\pi\right)^2\right]\\
 * \f]
 *
 * where \f$\sigma_+^{-2} = \beta\left|\cos\left(\frac{x_0}{4}\right)\right|\f$,
 * \f$\sigma_-^{-2} = \beta\left|\sin\left(\frac{x_0}{4}\right)\right|\f$ and
 * the normalisation constants
 *
 * \f[
 *    N_+(x_0) = \frac{1}{1+\rho(x_0)}\\
 *    N_-(x_0) = \frac{\rho(x_0)}{1+\rho(x_0)}
 * \f]
 *
 * The function \f$\rho\f$, which can take on values between 0 and 1 of
 * \f$x_0\in[0,\pi]\f$, is given by
 *
 * \f[
 *   \rho(x_0) = \left(\frac{\sigma_+}{\sigma_-}\right)^{3}
 * \exp\left[-4(\sigma_+^{-2}-\sigma_-^{-2})\right] \f]
 *
 * Note that it is sufficient to know \f$\tilde{p}(x|x_0,0)\f$ with
 * \f$x\in[-\pi,\pi]\f$ and \f$x_0\in[0,\pi]\f$
 *
 * \f$p(x|x_+,x_-)=p(x-x_-|x_+-x_-,0)\f$
 * \f$p(-x|-x_0,0) = p(x|x_0,0)\f$
 * \f$p(x|2\pi-x_0,0) = p(x|-x_0,0)\f$
 */

class ApproximateBesselProductDistribution {
public:
  /** @brief Constructor
   *
   * Create a new instance for a given value of \f$\beta\f$
   *
   * @param[in] beta_ Parameter \f$\beta\f$
   */
  ApproximateBesselProductDistribution(double beta_)
      : beta(beta_), distribution(0.0, 1.0), normal_distribution(0.0, 1.0),
        kmax(4) {}

  /** @brief Draw number from distribution for different \f$x_+,x_-\f$
   *
   * @param[in] engine Random number generator engine
   * @param[in] x_p Value of parameter \f$x_+\f$
   * @param[in] x_m Value of parameter \f$x_-\f$
   */
  template <class URNG>
  const double draw(URNG &engine, const double x_p, const double x_m) const {
    double x0 = x_p - x_m;
    // Flip sign if necessary to ensure dx is always non-negative
    // (this will be compensated for later)
    double sign_flip =
        (x0 < 0) ? -1 : +1; // Flip sign if difference is negative
    x0 *= sign_flip;
    if (x0 > M_PI) {
      x0 = 2. * M_PI - x0;
      sign_flip *= -1;
    }
    double N_p, sigma2_p_inv, sigma2_m_inv;
    compute_N_p_sigma2inv(beta, x0, N_p, sigma2_p_inv, sigma2_m_inv);
    double xi = distribution(engine);
    double sigma, xshift;
    if (xi <= N_p) {
      sigma = 1. / sqrt(sigma2_p_inv);
      xshift = 0.0;
    } else {
      sigma = 1. / sqrt(sigma2_m_inv);
      xshift = M_PI;
    }
    double x = sigma * normal_distribution(engine) + 0.5 * x0 - xshift;
    return mod_2pi(sign_flip * x + x_m);
  }

  /** @brief Evaluate distribution for a given value of \f$x\f$ and parameters
   * \f$x_+,x_-\f$
   *
   * Calculate the value of the distribution \f$p(x|x_+,x_-)\f$ at a point
   * \f$x\in[-\pi,\pi]\f$.
   *
   * @param[in] x Point \f$x\f$ at which to evaluate the distribution
   * @param[in] x_p Value of parameter \f$x_+\f$
   * @param[in] x_m Value of parameter \f$x_-\f$
   */
  double evaluate(const double x, const double x_p, const double x_m) const;

  /** @brief Return parameter \f$\beta\f$ */
  double get_beta() const { return beta; }

private:
  /** @brief Compute constants required for sampling from the distribution
   *
   * This method computes the quantities \f$N_+\f$ and \f$\sigma_\pm^{-2}\f$
   *
   * @param[in] beta Coupling constant \f$\beta\f$
   * @param[in] x0 Shift parameter \f$x_0\f$
   * @param[out] N_p Resulting probability of sampling from main peak
   * @param[out] sigma2_p_inv Resulting squared width \f$1/\sigma_+^2\f$ of
   * Gaussian
   * @param[out] sigma2_m_inv Resulting squared width \f$1/\sigma_-^2\f$ of
   * Gaussian
   */
  void compute_N_p_sigma2inv(const double beta, const double x0, double &N_p,
                             double &sigma2_p_inv, double &sigma2_m_inv) const;

  /** @brief parameter \f$\beta\f$*/
  const double beta;
  /** @brief Uniform distribution for sampling */
  mutable std::uniform_real_distribution<double> distribution;
  /** @brief Normal distribution for approximate sampling */
  mutable std::normal_distribution<double> normal_distribution;
  /** @brief Number  of terms in sum over periodic copies */
  const int kmax;
};

#endif // APPROXIMATEBESSELPRODUCTDISTRIBUTION
