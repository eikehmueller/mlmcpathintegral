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
