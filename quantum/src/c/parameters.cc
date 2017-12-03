/** @file parameters.cc
 * @brief Implementation of parameters.hh
 */
#include "parameters.hh"

void Parameters::readfile() {
  std::ifstream infile;
  std::string dummy;
  try {
    infile.open(filename,std::ios::in);
    if (not infile.is_open()) {
      std::cerr << "ERROR opening file \'" << filename << "\'" << std::endl;
      exit(-1);
    }
    // Read all parameters
    infile >> dummy >> M_lat;
    infile >> dummy >> T_final;
    infile >> dummy >> m0;
    infile >> dummy >> mu2;
    infile >> dummy >> perturbative;
    infile >> dummy >> n_burnin;
    infile >> dummy >> n_samples;
    infile.close();
  } catch (std::ifstream::failure e) {      
    std::cerr << "ERROR opening or reading file \'" << filename << "\': " << e.what() << std::endl;
  }
}

/** Print out all parameters */
void Parameters::show() {
  std::cout << " M_lat         = " << M_lat << std::endl;
  std::cout << " T_final       = " << T_final << std::endl;
  std::cout << " m0            = " << m0 << std::endl;
  std::cout << " mu2           = " << mu2 << std::endl;
  std::cout << " preturbative  = " << perturbative << std::endl;
  std::cout << " n_burnin      = " << n_burnin << std::endl;
  std::cout << " n_samples     = " << n_samples << std::endl;
}
