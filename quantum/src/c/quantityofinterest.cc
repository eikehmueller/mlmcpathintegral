#include "quantityofinterest.hh"
/** @file quantityofinterest.cc
 * @brief Implementation of quantityofinterest.hh
 */

/* Evaluate QoI */
const double QoIXsquared::evaluate(const Path* x_path) {
  double X2=0.0;
  for (unsigned int i=0; i<x_path->M_lat; ++i) {
    double tmp = x_path->data[i];
    X2 += tmp*tmp;
  }
  return X2/x_path->M_lat;
}

/* Evaluate QoI */
const double QoISusceptibility::evaluate(const Path* x_path) {
  double dx = x_path->data[0]-x_path->data[x_path->M_lat-1];
  double Q = mod_pi(dx);
  double chi = Q*Q;
  for (unsigned int i=1; i<x_path->M_lat; ++i) {
    dx = x_path->data[i]-x_path->data[i-1];
    Q = mod_pi(dx);
    chi += Q*Q;
  }
  return four_pi2_inv*chi/x_path->T_final;
}
