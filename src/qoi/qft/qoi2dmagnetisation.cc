#include "qoi2dmagnetisation.hh"
/** @file qoi2dmagnetisation.cc
 * @brief Implementation of qoi2dmagnetisation.hh
 */

/* Evaluate QoI */
const double QoI2DMagnetisation::evaluate(const std::shared_ptr<SampleState> phi_state) {
    if (phi_state->data.size() != 2*Mt_lat*Mx_lat) {
        mpi_parallel::cout << "ERROR: Evaluating QoI2DMagnetisation on state of wrong size." << std::endl;
        mpi_exit(EXIT_FAILURE);
    }
    // lambda function for working out linear index of link
    double mu[3] = {0.0, 0.0, 0.0};
    for (int i=0;i<Mt_lat;++i) {
        for (int j=0;j<Mx_lat;++j) {
            double theta = phi_state->data[2*lattice->vertex_cart2lin(i,j)];
            double phi = phi_state->data[2*lattice->vertex_cart2lin(i,j)+1];
            mu[0] += sin(theta)*cos(phi);
            mu[1] += sin(theta)*sin(phi);
            mu[2] += cos(theta);
        }
    }
    return (mu[0]*mu[0]+mu[1]*mu[1]+mu[2]*mu[2])/(Mt_lat*Mx_lat);
}
