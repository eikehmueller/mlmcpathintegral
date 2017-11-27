#include "quantityofinterest.hh"
/** @brief Implementation of quantityofinterest.hh
 */

/* Evaluate QoI */
const double QoIXsquared::evaluate(const Path* x_path) {
  return x_path->data[0]*x_path->data[0];
}
