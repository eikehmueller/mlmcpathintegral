#include "doublewellaction.hh"
/** @file doublewellaction.cc
 * @brief Implementation of doublewellaction.hh
 */

/** Evaluate action */
const double DoubleWellAction::evaluate(const Path* x_path) const {
  double ainv2 = 1./(a_lat*a_lat);
  double x_j = x_path->data[0];
  double x_j_squared = x_j*x_j;
  double x_diff = x_path->data[0]-x_path->data[M_lat-1];
  double S = 0.5*m0*(ainv2*x_diff*x_diff + mu2*x_j_squared) + lambda*exp(-0.5*inv_sigma2*x_j_squared);
  for (unsigned int j=1;j<M_lat;++j) {
    double x_j = x_path->data[j];
    double x_j_squared = x_j*x_j;
    double x_diff = x_j-x_path->data[j-1];
    S += 0.5*m0*(ainv2*x_diff*x_diff+mu2*x_j_squared)+lambda*exp(-0.5*inv_sigma2*x_j_squared);
  }
  return a_lat*S;
}

/** Calculate force */
void DoubleWellAction::force(const Path* x_path,
                             Path* p_path) const {
  double tmp_1 = m0/a_lat;
  double tmp_2 = 2.+a_lat*a_lat*mu2;
  double tmp_3 = a_lat*lambda*inv_sigma2;
  double X_j;
  X_j = x_path->data[0];
  p_path->data[0] = tmp_1*(tmp_2*X_j - x_path->data[M_lat-1] - x_path->data[1])-tmp_3*X_j*exp(-0.5*inv_sigma2*X_j*X_j);
  // Interior points
  for (unsigned int j=1;j<M_lat-1;++j) {
    X_j = x_path->data[j];
    p_path->data[j] = tmp_1*(tmp_2*X_j - x_path->data[j-1] - x_path->data[j+1])-tmp_3*X_j*exp(-0.5*inv_sigma2*X_j*X_j);
  }
  X_j = x_path->data[M_lat-1];
  p_path->data[M_lat-1] = tmp_1*(tmp_2*X_j - x_path->data[M_lat-2] - x_path->data[0])-tmp_3*X_j*exp(-0.5*inv_sigma2*X_j*X_j);
}
