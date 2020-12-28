#include "qoiavgplaquette.hh"
/** @file qoiavgplaquette.cc
 * @brief Implementation of qoiavgplaquette.hh
 */

/* Evaluate QoI */
const double QoIAvgPlaquette::evaluate(const std::shared_ptr<SampleState> phi_state) {
    // lambda function for working out linear index of link
    double S_plaq = 0.0;
    for (int i=0;i<Mt_lat;++i) {
        for (int j=0;j<Mx_lat;++j) {
            double theta = phi_state->data[lattice->link_cart2lin(i,  j,  0)]
                         + phi_state->data[lattice->link_cart2lin(i+1,j,  1)]
                         - phi_state->data[lattice->link_cart2lin(i,  j+1,0)]
                         - phi_state->data[lattice->link_cart2lin(i,  j  ,1)];
            S_plaq += cos(theta);
        }
    }
    return S_plaq/(Mx_lat*Mt_lat);
}
