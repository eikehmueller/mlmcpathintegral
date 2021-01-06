#include "qoixsquared.hh"
/** @file qoixsquared.cc
 * @brief Implementation of qoixsquared.hh
 */

/* Evaluate QoI */
const double QoIXsquared::evaluate(const std::shared_ptr<SampleState> x_path) {
    double X2=0.0;
    for (unsigned int i=0; i<M_lat; ++i) {
        double tmp = x_path->data[i];
        X2 += tmp*tmp;
    }
    return X2/M_lat;
}