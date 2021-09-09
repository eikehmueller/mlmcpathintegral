#include "qftaction.hh"
/** @file qftaction.cc
 * @brief Implementation of qftaction.hh
 */

/* Return lattice information */
std::string QFTAction::info_string() const {
    std::stringstream sstr;
    sstr << "lattice = " << lattice->getMt_lat() << " x " << lattice->getMx_lat();
    return sstr.str();
}
