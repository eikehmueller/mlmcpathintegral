#include "path.hh"
/** @file path.cc
 * @brief Implementation of path.hh
 */

/* Copy data */
void Path::copy(const std::shared_ptr<Path> other) {
  assert(other->M_lat == M_lat);
  std::copy_n(other->data,M_lat,data);
}

/* Copy data with stride */
void Path::copy_strided(const std::shared_ptr<Path> other,
                        const int stride) {
  unsigned int s = abs(stride);
  if (stride > 0) {
    assert(other->M_lat == s*M_lat);
    for (unsigned int j=0;j<M_lat;++j) {
      data[j] = other->data[s*j];
    }
  } else {
    assert(s*other->M_lat == M_lat);
    for (unsigned int j=0;j<M_lat/s;++j) {
      data[s*j] = other->data[j];
    }
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
