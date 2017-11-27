#include "quantityofinterest.hh"
/** @brief Implementation of quantityofinterest.hh
 */

/* Evaluate QoI */
const double QoIXsquared::evaluate(const double* x_path) {
  return x_path[0]*x_path[0];
}
