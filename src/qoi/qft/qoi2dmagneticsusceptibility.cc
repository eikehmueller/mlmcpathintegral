#include "qoi2dmagneticsusceptibility.hh"
/** @file qoi2dmagneticsusceptibility.cc
 * @brief Implementation of qoi2ddmagneticsusceptibility.hh
 */

/* Evaluate QoI */
const double QoI2DMagneticSusceptibility::evaluate(const std::shared_ptr<SampleState> phi_state) {
    // lambda function for working out linear index of link
    double mu[3] = {0.0, 0.0, 0.0};
    for (unsigned int ell=0;ell<lattice->getNvertices();++ell) {
        double theta = phi_state->data[2*ell];
        double phi = phi_state->data[2*ell+1];
        double sin_theta = sin(theta);
        mu[0] += sin_theta*cos(phi);
        mu[1] += sin_theta*sin(phi);
        mu[2] += cos(theta);                
    }
    return (mu[0]*mu[0]+mu[1]*mu[1]+mu[2]*mu[2])/lattice->getNvertices();
}
