#include "renormalisation.hh"

/** @file renormalisation.cc
 * @brief Implementation of renormalisation.hh 
 */

double  RenormalisedRotorParameters::deltaI(const double xi) {
  double a = 1.-17./8.*xi*Sigma_hat(xi,2)+1./2.*xi*Sigma_hat(xi,4);
  double b = 1.-xi*Sigma_hat(xi,2)*(2.+xi*Sigma_hat(xi,2))+xi*xi*Sigma_hat(xi,4);
  return 0.5*a/b;
}

double  RenormalisedRotorParameters::Sigma_hat(const double xi,
                                               const unsigned int p) {
  if (p % 2 == 0) {
    if (p == 0) {
      // Ratio is one if p==0
      return 1.0;
    } else {
      unsigned int mmax = 100; // Number of terms in sum
      double num = 0.0;
      double denom = 1.0;
      for (unsigned int m=1;m<mmax;++m) {
        double exp_factor = exp(-0.5*xi*m*m);
        num += 2.*pow(m,p)*exp_factor;
        denom += 2.*exp_factor;
      }
      return num/denom;
    }
  } else {
    // Sum is zero if p is odd (symmetry)
    return 0.0;
  }
}
