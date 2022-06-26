#include "qoisusceptibility.hh"
/** @file qoisusceptibility.cc
 * @brief Implementation of quantityofinterest.hh
 */

/* Evaluate QoI */
const double
QoISusceptibility::evaluate(const std::shared_ptr<SampleState> x_path) {
  if (x_path->data.size() != M_lat) {
    mpi_parallel::cout
        << "ERROR: Evaluating QoISusceptibility on path of wrong size."
        << std::endl;
    mpi_exit(EXIT_FAILURE);
  }
  double dx = x_path->data[0] - x_path->data[M_lat - 1];
  double Q = mod_2pi(dx);
  for (unsigned int i = 1; i < M_lat; ++i) {
    dx = x_path->data[i] - x_path->data[i - 1];
    Q += mod_2pi(dx);
  }
  double chi = Q * Q;
  return four_pi2_inv * chi / T_final;
}
