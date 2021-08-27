#include "qoi2dmagneticsusceptibility.hh"
/** @file qoi2dmagneticsusceptibility.cc
 * @brief Implementation of qoi2ddmagneticsusceptibility.hh
 */

/* Evaluate QoI */
const double QoI2DMagneticSusceptibility::evaluate(const std::shared_ptr<SampleState> phi_state) {
    // lambda function for working out linear index of link
    double mu[3] = {0.0, 0.0, 0.0};
    for (int i=0;i<Mt_lat;++i) {
        for (int j=0;j<Mx_lat;++j) {
            if ( (not rotated) or ((i+j)%2==0) ) {
                double theta;
                double phi;
                if (rotated) {
                    theta = phi_state->data[2*lattice->diag_vertex_cart2lin(i,j)];
                    phi = phi_state->data[2*lattice->diag_vertex_cart2lin(i,j)+1];                    
                } else {
                    theta = phi_state->data[2*lattice->vertex_cart2lin(i,j)];
                    phi = phi_state->data[2*lattice->vertex_cart2lin(i,j)+1];                    
                }
                mu[0] += sin(theta)*cos(phi);
                mu[1] += sin(theta)*sin(phi);
                mu[2] += cos(theta);                
            }
        }
    }
    int N_lat;
    if (rotated) {
        N_lat = Mt_lat*Mx_lat/2;
    } else {
        N_lat = Mt_lat*Mx_lat;
    }
    return (mu[0]*mu[0]+mu[1]*mu[1]+mu[2]*mu[2])/N_lat;
}
