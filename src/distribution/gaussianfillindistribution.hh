#ifndef GAUSSIANFILLINDISTRIBUTION_HH
#define GAUSSIANFILLINDISTRIBUTION_HH GAUSSIANFILLINDISTRIBUTIONHH
#include <random>
#include <set>
#include <array>
#include <vector>
#include <algorithm>
#include <cmath>
#include <iostream>
#include "common/auxilliary.hh"

/** @file gaussianfillindistribution.hh
 *
 * @brief Header file for preconditioned 4D Gaussian distribution used for fill-in in 
 *        quenched Schwinger model
 */

/** @class GaussianFillinDistribution
 *
 * @brief Class for sampling from a Gaussian approximation of the fill-in distribution
 *
 * During fill-in for the quenched Schwinger model we need to sample from the following
 * distribution:
 * 
 * \f[
 * p_0(\theta|\phi) \propto \exp\left[ \beta \left(  \cos(\theta_1-\theta_2-\phi_{12})
 *                                                 + \cos(\theta_2-\theta_3-\phi_{23})
 *                                                 + \cos(\theta_3-\theta_4-\phi_{34})
 *                                                 + \cos(\theta_4-\theta_1-\phi_{41})
 *                               \right) \right]
 * \f]
 * for given values of \f$\phi\f$. This class approximates this distribution by a 
 * multivariate Gaussian.
 * 
 * If the parameter add_gaussian_noise is set to false, the two Gaussian distributions
 * are further replaced by delta-functions, in other words only one of the most likely
 * configurations is returned.
 */

class GaussianFillinDistribution {
public:
  /** @brief Constructor
   * 
   * Create a new instance for a given value of \f$\beta\f$
   *
   * @param[in] beta_ Parameter \f$\beta\f$
   * @param[in] add_gaussian_noise_ Add Gaussian noise? If false, only take peak value
   */
  GaussianFillinDistribution(const double beta_,
                             const bool add_gaussian_noise_ = false) :
    beta(beta_),
    add_gaussian_noise(add_gaussian_noise_),
    uniform_distribution(0.0,1.0),
    normal_distribution(0.0,1.0),
    sqrt2(sqrt(2.0)) {
        if (not add_gaussian_noise) {
            mpi_parallel::cerr << "ERROR: sampling only from peak currently broken." << std::endl;
            mpi_exit(EXIT_FAILURE);
        }
        int n_offsets = 1;
        // For large values of beta the peaks in period copies of the unit cell
        // are strongly suppressed.
        if (beta > 72.0) n_offsets = 0;
        construct_peaks(n_offsets);
    }

  /** @brief Draw number directly from Gaussian distribution for given \f$\phi\f$
   *
   * @param[in] engine Random number generator engine
   * @param[in] phi_12 Parameter \f$\phi_{12}\f$
   * @param[in] phi_23 Parameter \f$\phi_{23}\f$
   * @param[in] phi_34 Parameter \f$\phi_{34}\f$
   * @param[in] phi_41 Parameter \f$\phi_{41}\f$
   * @param[out] theta_1 Resulting value \f$\theta_1\f$
   * @param[out] theta_2 Resulting value \f$\theta_2\f$
   * @param[out] theta_3 Resulting value \f$\theta_3\f$
   * @param[out] theta_4 Resulting value \f$\theta_4\f$
   */
  template <class URNG>
  void draw(URNG& engine,
            const double phi_12, const double phi_23,
            const double phi_34, const double phi_41,
            double& theta_1, double& theta_2,
            double& theta_3, double& theta_4) {
      double eta_1, eta_2, eta_3; // Coordinates in non-trivial 3d subspace
      double sigma; // Variance of Gaussians       
      double Phi = 0.25*(phi_12+phi_23+phi_34+phi_41);
      double Phi_star = Phi;
      // Map Phi to the range [0,pi/2]
      bool swap_eta = false;
      bool shift_eta = false;
      if (Phi_star < 0) {
          Phi_star *= -1.0;
          swap_eta = true;
      }
      if (Phi_star > 0.5*M_PI) {
          Phi_star = M_PI - Phi_star;
          swap_eta = not(swap_eta);
          shift_eta = true;
      }
      // Probability of sampling from main Gaussian
      double p_c = get_pc(Phi_star);
      double xi = uniform_distribution(engine);
      if (xi < p_c) {
          // Draw from main peak
          eta_1 = 0.0;
          eta_2 = 0.0;
          eta_3 = 0.0;
          sigma = 1./sqrt(4.*beta*cos(Phi_star));
      } else {
          // Draw from secondary peak           
          eta_1 = M_PI;
          eta_2 = 0.0;
          eta_3 = 0.5*M_PI;
          sigma = 1./sqrt(4.*beta*sin(Phi_star));
      }
      if (add_gaussian_noise) {
          eta_1 += sqrt2*sigma*normal_distribution(engine);
          eta_2 += sqrt2*sigma*normal_distribution(engine);
          eta_3 += sigma*normal_distribution(engine);
      }
      if (swap_eta) {
          std::swap(eta_1,eta_2);
      }
      if (shift_eta) {
          eta_1 += M_PI;
          eta_2 += M_PI;
      }
      // Draw uniform shift
      double omega = 2.*M_PI*uniform_distribution(engine);
      theta_1 = mod_2pi(0.5*(+eta_1+eta_2+eta_3) + omega);
      theta_2 = mod_2pi(0.5*(+eta_1-eta_2-eta_3) + omega + Phi-phi_12);
      theta_3 = mod_2pi(0.5*(-eta_1-eta_2+eta_3) + omega + 2.*Phi-phi_12-phi_23);
      theta_4 = mod_2pi(0.5*(-eta_1+eta_2-eta_3) + omega + 3.*Phi-phi_12-phi_23-phi_34);
    }
  
  /** @brief Evaluate distribution for a given value of \f$\theta\f$
   *         and parameters \f$\phi\f$
   *
   * Calculate the value of the distribution \f$p(\theta|\phi)\f$ at a point
   * \f$\theta\f$.
   *
   * @param[in] theta_1 Value of \f$\theta_1\f$
   * @param[in] theta_2 Value of \f$\theta_2\f$
   * @param[in] theta_3 Value of \f$\theta_3\f$
   * @param[in] theta_4 Value of \f$\theta_4\f$
   * @param[in] phi_12 Value of \f$\phi_{12}\f$
   * @param[in] phi_23 Value of \f$\phi_{23}\f$
   * @param[in] phi_34 Value of \f$\phi_{34}\f$
   * @param[in] phi_41 Value of \f$\phi_{41}\f$
   */
  double evaluate(const double theta_1, const double theta_2,
                  const double theta_3, const double theta_4,
                  const double phi_12, const double phi_23,
                  const double phi_34, const double phi_41) const;
  
    /** @brief Return parameter \f$\beta\f$ */
  double get_beta() const { return beta; }
private:

  /** @brief Construct peaks used the density in the Gaussian approximation
   *
   * @param[in] n_offsets Number of offset copies of unit cell
   */
  void construct_peaks(int n_offsets);


  /** @brief Compute probability of sampling from main Gaussian
   *
   * @param[in] Phi angle \f$\Phi\f$
   */
  double get_pc(const double Phi) const {
    if (Phi < 0.125*M_PI) {
        return 1.0;
    } else if (Phi > 0.375*M_PI) {
        return 0.0;
    } else {
        double sigma2_p_inv = beta*cos(Phi);
        double sigma2_m_inv = beta*sin(Phi);
        double rho = pow(sigma2_p_inv/sigma2_m_inv,1.5)*exp(-4.0*(sigma2_p_inv-sigma2_m_inv));
        return 1./(1.+rho);        
    }
  }
    
  /** @brief parameter \f$\beta\f$*/
  const double beta;
  /** @brief Uniform distribution for shift */
  mutable std::uniform_real_distribution<double> uniform_distribution;
  /** @brief Normal distribution for approximate sampling */
  mutable std::normal_distribution<double> normal_distribution;
  /** @brief Precomputed value of \f$\sqrt{2}\f$ */
  const double sqrt2;
  /** @brief Add Gaussian noise? If false, only take peak value */
  const bool add_gaussian_noise;
  /** @brief Locations of main peaks in box and nearest neighbours */
  mutable std::vector<std::array<double,3> > main_peaks;
  /** @brief Locations of secondary peaks in box and nearest neighbours */
  mutable std::vector<std::array<double,3> > secondary_peaks;
};

#endif // GAUSSIANFILLINDISTRIBUTION
