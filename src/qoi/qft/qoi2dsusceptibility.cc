#include "qoi2dsusceptibility.hh"
/** @file qoi2dsusceptibility.cc
 * @brief Implementation of qoi2dsusceptibility.hh
 */

/* Evaluate QoI */
const double QoI2DSusceptibility::evaluate(const std::shared_ptr<SampleState> phi_state) {
    // lambda function for working out linear index of link
    double Q = 0.0;
    for (int i=0;i<Mt_lat;++i) {
        for (int j=0;j<Mx_lat;++j) {
            double theta = phi_state->data[lattice->link_cart2lin(i,  j,  0)]
                         + phi_state->data[lattice->link_cart2lin(i+1,j,  1)]
                         - phi_state->data[lattice->link_cart2lin(i,  j+1,0)]
                         - phi_state->data[lattice->link_cart2lin(i,  j  ,1)];
            Q += mod_2pi(theta);
        }
    }
    return four_pi2_inv*Q*Q;
}
