#include "path.hh"
/** @file path.cc
 * @brief Implementation of path.hh
 */

/* Copy data */
void Path::copy(const std::shared_ptr<Path> other) {
  assert(other->M_lat == M_lat);
  std::copy_n(other->data,M_lat,data);
}

/* Copy coarse data from path on coarser level */
void Path::copy_from_coarse(const std::shared_ptr<Path> other) {
  assert(2*other->M_lat == M_lat);
  for (unsigned int j=0;j<M_lat/2;++j) {
    data[2*j] = other->data[j];
  }
}

/* Copy coarse data from path on finer level */
void Path::copy_from_fine(const std::shared_ptr<Path> other) {
  assert(other->M_lat == 2*M_lat);
  for (unsigned int j=0;j<M_lat;++j) {
    data[j] = other->data[2*j];
  }
}


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
