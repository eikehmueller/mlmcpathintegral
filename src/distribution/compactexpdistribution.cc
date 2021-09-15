#include "compactexpdistribution.hh"
/** @file compactexpdistribution.cc
 * @brief Implementation of compactexpdistribution.hh
 */

/** Evaluate distribution */
double CompactExpDistribution::evaluate(const double x, const double sigma) const {
    return 0.5*sigma*exp(sigma*x)/sinh(sigma);
}
