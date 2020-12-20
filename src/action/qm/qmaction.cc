#include "qmaction.hh"
/** @file qmaction.cc
 * @brief Implementation of qmaction.hh
 */

/* Copy coarse data from path on coarser level */
void QMAction::copy_from_coarse(const std::shared_ptr<SampleState> x_coarse,
                                std::shared_ptr<SampleState> x_path) {
    assert(M_lat == x_path->data.size());
    assert(M_lat == 2*x_coarse->data.size());
    for (unsigned int j=0; j<M_lat/2; ++j) {
        x_path->data[2*j] = x_coarse->data[j];
    }
}

/* Copy coarse data from path on finer level */
void QMAction::copy_from_fine(const std::shared_ptr<SampleState> x_fine,
                              std::shared_ptr<SampleState> x_path) {
    assert(M_lat == x_path->data.size());
    assert(2*M_lat == x_fine->data.size());
    for (unsigned int j=0; j<M_lat; ++j) {
        x_path->data[j] = x_fine->data[2*j];
    }
}

void QMAction::check_coarsening_is_permitted(const unsigned int n_level) {
    if ( (M_lat>>n_level)<<n_level == M_lat) {
        mpi_parallel::cout << "M_lat = " << M_lat << " = 2^{" << n_level << "-1} * " << (M_lat>>(n_level-1)) << std::endl;
    } else {
        mpi_parallel::cout << "ERROR: M_lat = " << M_lat << " is not a multiple of 2^{n_level} = 2^{"<<n_level << "}" << std::endl;
    }
}
