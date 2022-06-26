#include "samplestate.hh"
/** @file samplestate.cc
 * @brief Implementation of samplestate.hh
 */

/** Save state to disk */
void SampleState::save_to_disk(const std::string filename) {
  std::ofstream file;
  file.open(filename.c_str());
  file << M << std::endl;
  for (unsigned int i = 0; i < M; ++i) {
    file << data[i] << " ";
  }
  file << std::endl;
  file.close();
}
