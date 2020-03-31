#include "renormalisation.hh"

/** @file renormalisation.cc
 * @brief Implementation of renormalisation.hh 
 */

double  RenormalisedRotorParameters::deltaI(const double xi) {
  double xi = T_final/m0;
  double S_hat2 = Sigma_hat(xi,2);
  double S_hat4 = Sigma_hat(xi,4);
  return 0.5*(1.-2.*xi*S_hat2+0.5*xi*xi*(S_hat4-S_hat2*S_hat2))/(1.-2.*xi*S_hat2+xi*xi*(S_hat4-S_hat2*S_hat2));
}
