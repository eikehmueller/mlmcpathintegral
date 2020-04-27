#include "path.hh"
/** @file path.cc
 * @brief Implementation of path.hh
 */

/** Save path to disk */
void Path::save_to_disk(const std::string filename) {
  std::ofstream file;
  file.open(filename.c_str());
  file << M_lat << " " << T_final << std::endl;
  for (unsigned int i=0; i<M_lat; ++i) {
    file << data[i] << " ";
  }
  file << std::endl;
  file.close();
}
