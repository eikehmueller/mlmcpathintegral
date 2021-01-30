#include "auxilliary.hh"
/** @file auxilliary.cc
 * @brief Implementation of auxilliary.hh
 */

/* Implementation of Sigma_hat */
double Sigma_hat(const double xi,
                 const unsigned int p) {
    if (p % 2 == 0) {
        if (p == 0) {
            // Ratio is one if p==0
            return 1.0;
        } else {
            unsigned int mmax = 100; // Number of terms in sum
            double num = 0.0;
            double denom = 1.0;
            for (unsigned int m=1; m<mmax; ++m) {
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

/* Implementation of log_factorial */
double log_factorial(unsigned int n) {
  double s=0.0;
  for (unsigned int k=2;k<=n;++k) {
    s += log(k);
  }
  return s;
}

/* Implementation of log_nCk */
double log_nCk(unsigned int n, unsigned int k) {
 return log_factorial(n)-log_factorial(k)-log_factorial(n-k);
}
