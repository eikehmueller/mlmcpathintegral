#include "qftaction.hh"
/** @file qftaction.cc
 * @brief Implementation of qftaction.hh
 */

/* Check whether coarsening is permitted */
void QFTAction::check_coarsening_is_permitted(const unsigned int n_level) {
    if ( (Mt_lat>>(n_level-1))<<(n_level-1) == Mt_lat) {
        mpi_parallel::cout << "M_{t,lat} = " << Mt_lat << " = 2^{" << n_level << "-1} * " << (Mt_lat>>(n_level-1)) << std::endl;
    } else {
        mpi_parallel::cout << "ERROR: M_{t,lat} = " << Mt_lat << " is not a multiple of 2^{n_level-1} = 2^{"<<n_level-1 << "}" << std::endl;
        mpi_exit(-1);
    }
    if ( (Mx_lat>>(n_level-1))<<(n_level-1) == Mx_lat) {
        mpi_parallel::cout << "M_{x,lat} = " << Mx_lat << " = 2^{" << n_level << "-1} * " << (Mx_lat>>(n_level-1)) << std::endl;
    } else {
        mpi_parallel::cout << "ERROR: M_{x,lat} = " << Mx_lat << " is not a multiple of 2^{n_level-1} = 2^{"<<n_level-1 << "}" << std::endl;
        mpi_exit(-1);
    }
}
