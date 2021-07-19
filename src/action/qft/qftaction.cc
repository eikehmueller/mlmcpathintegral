#include "qftaction.hh"
/** @file qftaction.cc
 * @brief Implementation of qftaction.hh
 */

/* Return lattice information */
std::string QFTAction::info_string() const {
    std::stringstream sstr;
    sstr << "lattice = " << Mt_lat << " x " << Mx_lat;
    return sstr.str();
}
