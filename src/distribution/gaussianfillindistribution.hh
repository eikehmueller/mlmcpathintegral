#ifndef GAUSSIANFILLINDISTRIBUTION_HH
#define GAUSSIANFILLINDISTRIBUTION_HH GAUSSIANFILLINDISTRIBUTIONHH
#include <random>
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
 */

class GaussianFillinDistribution {
public:
  /** @brief Constructor
   * 
   * Create a new instance for a given value of \f$\beta\f$
   *
   * @param[in] beta_ Parameter \f$\beta\f$
   * @param[in] alpha_pcn_ pCN offset parameter \f$\alpha\f$
   */
  GaussianFillinDistribution(const double beta_, const double alpha_pcn_ = 0.9) :
    beta(beta_),
    alpha_pcn(alpha_pcn_),
    alpha_pcn_comp(sqrt(1.-alpha_pcn_*alpha_pcn_)),
    uniform_distribution(0.0,1.0),
    normal_distribution(0.0,1.0),
    sqrt2(sqrt(2.0)) {}

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
      double mu_1, mu_2, mu_3, sigma; // Mean and variance of Gaussians       
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
      /* === DEBUG === */
      // Only pick a single peak (to be consistent with pCN proposal)
      if (p_c > 0.5) {
          p_c = 1.0;
      } else {
          p_c = 0.0;
      }
      /* === END DEBUG === */
      double xi = uniform_distribution(engine);
      if (xi < p_c) {
          // Draw from main peak
          mu_1 = 0.0;
          mu_2 = 0.0;
          mu_3 = 0.0;
          sigma = 1./sqrt(4.*beta*cos(Phi_star));
      } else {
          // Draw from secondary peak           
          mu_1 = M_PI;
          mu_2 = 0.0;
          mu_3 = 0.5*M_PI;
          sigma = 1./sqrt(4.*beta*sin(Phi_star));
      }
      double eta_1 = mu_1 + sqrt2*sigma*normal_distribution(engine);
      double eta_2 = mu_2 + sqrt2*sigma*normal_distribution(engine);
      double eta_3 = mu_3 + sigma*normal_distribution(engine);
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
    
  /** @brief Draw using preconditioned Crank Nicolson method
   *
   * @param[in] engine Random number generator engine
   * @param[in] phi_n_12 Parameter \f$\phi_{12}\f$ at previous step
   * @param[in] phi_n_23 Parameter \f$\phi_{23}\f$ at previous step
   * @param[in] phi_n_34 Parameter \f$\phi_{34}\f$ at previous step
   * @param[in] phi_n_41 Parameter \f$\phi_{41}\f$ at previous step
   * @param[in] theta_n_1 Value \f$\theta_1\f$ at previous step
   * @param[in] theta_n_2 Value \f$\theta_2\f$ at previous step
   * @param[in] theta_n_3 Value \f$\theta_3\f$ at previous step
   * @param[in] theta_n_4 Value \f$\theta_4\f$ at previous step
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
  void draw_pcn(URNG& engine,
                const double phi_n_12, const double phi_n_23,
                const double phi_n_34, const double phi_n_41,
                const double theta_n_1, const double theta_n_2,
                const double theta_n_3, const double theta_n_4,
                const double phi_12, const double phi_23,
                const double phi_34, const double phi_41,
                double& theta_1, double& theta_2,
                double& theta_3, double& theta_4) {
      // === Step 1 ===
      // Construct three dimensional eta_n for previous step
      double eta_1 = 0.5*(theta_n_1+theta_n_2-theta_n_3-theta_n_4) + 0.5*(phi_n_41-phi_n_23);
      double eta_2 = 0.5*(theta_n_1-theta_n_2-theta_n_3+theta_n_4) + 0.5*(phi_n_34-phi_n_12);
      double eta_3 = 0.5*(theta_n_1-theta_n_2+theta_n_3-theta_n_4) + 0.25*(-phi_n_12+phi_n_23-phi_n_34+phi_n_41);
      double Phi_n = 0.25*(phi_n_12+phi_n_23+phi_n_34+phi_n_41);
      // === Step 2 ===
      // Map Phi_n into range [0,pi/2]
      bool swap_eta = false;
      double Phi_n_star = Phi_n;
      if (Phi_n_star < 0.) {
          Phi_n_star *= -1.0;
          swap_eta = true;
      }
      if (Phi_n_star > 0.5*M_PI) {
          Phi_n_star = M_PI-Phi_n_star;
          swap_eta = not(swap_eta);
          eta_1 = eta_1+M_PI;
          eta_2 = eta_2+M_PI;
      }
      if (swap_eta) {
          std::swap(eta_1,eta_2);        
      }
      // === Step 3 ===
      // Rescale eta_n with covariance from previous step to map it to
      // variable drawn from a 3d normal N(0,Id)
      double p_c_n = get_pc(Phi_n_star);
      double d_eta_2, d_eta_3;
      double sigma_inv_n;
      if (p_c_n > 0.5) {
          // Main peak
          d_eta_2 = 0.0;
          d_eta_3 = 0.0;
          sigma_inv_n = sqrt(2.*beta*cos(Phi_n_star));
      } else {
          // Secondary peak
          d_eta_2 = M_PI;
          d_eta_3 = 0.5*M_PI;
          sigma_inv_n = sqrt(2.*beta*sin(Phi_n_star));         
      }
      // Account for symmetry under shifts by integer multiples of 
      // (\pi,\pi,\pi), (\pi,-\pi,\pi) and (\pi,-\pi,-\pi)
      double eta_1_bar = mod_pi(0.5*(eta_1+eta_2+d_eta_2));
      double eta_2_bar = mod_pi(0.5*(eta_3+d_eta_3-(eta_2+d_eta_2)));
      double eta_3_bar = mod_pi(0.5*(eta_1-(eta_3+d_eta_3)));
      eta_1 = eta_1_bar + eta_2_bar + eta_3_bar;
      eta_2 = eta_1_bar - eta_2_bar - eta_3_bar;
      eta_3 = eta_1_bar + eta_2_bar - eta_3_bar;      
      eta_1 = sigma_inv_n*mod_2pi(eta_1);
      eta_2 = sigma_inv_n*mod_2pi(eta_2);
      eta_3 = sqrt2*sigma_inv_n*mod_2pi(eta_3);
      // === Step 4 ===
      // pCN draw
      eta_1 = alpha_pcn_comp*eta_1 + alpha_pcn*normal_distribution(engine);
      eta_2 = alpha_pcn_comp*eta_2 + alpha_pcn*normal_distribution(engine);
      eta_3 = alpha_pcn_comp*eta_3 + alpha_pcn*normal_distribution(engine);
      // === Step 5 ===
      // Map Phi to range [0,pi/2] and rescale eta_n with covariance from current step
      double Phi = 0.25*(phi_12+phi_23+phi_34+phi_41);
      double Phi_star = Phi;
      bool shift_eta = false;
      if (Phi_star < 0.) {
          Phi_star *= -1.0;
          swap_eta = true;
      } else {
          swap_eta = false;
      }
      if (Phi_star > 0.5*M_PI) {
          Phi_star = M_PI-Phi_star;
          swap_eta = not(swap_eta);          
          shift_eta = true;          
      }
      double sigma;
      double p_c = get_pc(Phi_star);
      if (p_c > 0.5) {
          d_eta_2 = 0.0;
          d_eta_3 = 0.0;
          sigma = 1./sqrt(2.*beta*cos(Phi_star));
      } else {
          d_eta_2 = M_PI;
          d_eta_3 = 0.5*M_PI;
          sigma = 1./sqrt(2.*beta*sin(Phi_star));         
      }
      eta_1 = sigma*eta_1;
      eta_2 = sigma*eta_2-d_eta_2;
      eta_3 = sigma/sqrt2*eta_3-d_eta_3;
      if (swap_eta) {
          std::swap(eta_1,eta_2);        
      }
      if (shift_eta) {
        eta_1 = eta_1+M_PI;
        eta_2 = eta_2+M_PI;
      }
      // === Step 6 ===
      // Draw uniform random variable
      double omega = 2.*M_PI*uniform_distribution(engine);
      // === Step 7 ===
      // Map back to the original 4d space
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

  /** @brief Compute probability of sampling from main Gaussian
   *
   * @param[in] Phi angle \f$\Phi\f$
   */
  double get_pc(const double Phi) const {
    double sigma2_p_inv = beta*cos(Phi);
    double sigma2_m_inv = beta*sin(Phi);
    double rho = pow(sigma2_m_inv/sigma2_p_inv,1.5)*exp(-4.0*(sigma2_p_inv-sigma2_m_inv));
    return 1./(1.+rho);
  }
    
  /** @brief parameter \f$\beta\f$*/
  const double beta;
  /** @brief Uniform distribution for shift */
  mutable std::uniform_real_distribution<double> uniform_distribution;
  /** @brief Normal distribution for approximate sampling */
  mutable std::normal_distribution<double> normal_distribution;
  /** @brief Precomputed value of \f$\sqrt{2}\f$ */
  const double sqrt2;
  /** @brief pCN parameter \f$\alpha\f$ */
  const double alpha_pcn;
  /** @brief complementary pCN parameter \f$\sqrt{1-\alpha^2}\f$ */
  const double alpha_pcn_comp;
};

#endif // GAUSSIANFILLINDISTRIBUTION
