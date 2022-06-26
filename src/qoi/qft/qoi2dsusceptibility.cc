#include "qoi2dsusceptibility.hh"
/** @file qoi2dsusceptibility.cc
 * @brief Implementation of qoi2dsusceptibility.hh
 */

/* Evaluate QoI */
const double
QoI2DSusceptibility::evaluate(const std::shared_ptr<SampleState> phi_state) {
  if (phi_state->data.size() != 2 * Mt_lat * Mx_lat) {
    mpi_parallel::cout
        << "ERROR: Evaluating QoI2DSusceptibility on state of wrong size."
        << std::endl;
    mpi_exit(EXIT_FAILURE);
  }
  // lambda function for working out linear index of link
  double Q = 0.0;
  for (int i = 0; i < Mt_lat; ++i) {
    for (int j = 0; j < Mx_lat; ++j) {
      double theta = phi_state->data[lattice->link_cart2lin(i, j, 0)] +
                     phi_state->data[lattice->link_cart2lin(i + 1, j, 1)] -
                     phi_state->data[lattice->link_cart2lin(i, j + 1, 0)] -
                     phi_state->data[lattice->link_cart2lin(i, j, 1)];
      Q += mod_2pi(theta);
    }
  }
  return four_pi2_inv * Q * Q;
}

/* Analytical expression for topological susceptibility */
double quenchedschwinger_chit_analytical(const double beta,
                                         const unsigned int n_plaq) {
  return n_plaq / beta * Phi_chit(beta, n_plaq);
}

/* Perturbative expression for topological susceptibility */
double quenchedschwinger_chit_perturbative(const double beta,
                                           const unsigned int n_plaq) {
  return n_plaq / beta * Phi_chit_perturbative(beta, n_plaq);
}

/* Analytical expression for variance of topological susceptibility in the
 * continuum limit */
double
quenchedschwinger_var_chit_continuum_analytical(const double beta,
                                                const unsigned int n_plaq) {
  double zeta = 4 * M_PI * M_PI * beta / n_plaq;
  double Sigma_hat_2 = Sigma_hat(zeta, 2);
  double Sigma_hat_4 = Sigma_hat(zeta, 4);
  return Sigma_hat_4 - Sigma_hat_2 * Sigma_hat_2;
}
