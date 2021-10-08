#include "qoi2dphisquared.hh"
/** @file qoi2dphisquared.cc
 * @brief Implementation of qoi2dphisquared.hh
 */

/* Evaluate QoI */
const double QoI2DPhiSquared::evaluate(const std::shared_ptr<SampleState> phi_state) {
    double phi_squared = 0.0;
    for (int ell=0;ell<M_lat;++ell) {
        double phi_n = phi_state->data[ell];
        phi_squared += phi_n*phi_n;
    }
    return phi_squared/M_lat;
}