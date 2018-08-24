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
    int renormalisation_int;
    int sampler_int;
    // Read all parameters
    infile >> dummy >> M_lat;
    infile >> dummy >> T_final;
    infile >> dummy >> m0;
    infile >> dummy >> mu2;
    infile >> dummy >> lambda;
    infile >> dummy >> sigma;
    infile >> dummy >> renormalisation_int;
    infile >> dummy >> n_burnin;
    infile >> dummy >> n_samples;
    infile >> dummy >> sampler_int;
    infile >> dummy >> T_hmc;
    infile >> dummy >> dt_hmc;
    infile >> dummy >> n_burnin_hmc;
    infile.close();
    switch (renormalisation_int) {
      case 0: renormalisation = RenormalisationNone; break;
      case 1: renormalisation = RenormalisationPerturbative; break;
      case 2: renormalisation = RenormalisationExact; break;
      default: {
        std::cout << " ERROR: Unknown renormalisation: " << renormalisation_int << std::endl;
        std::cout << "        allowed values are 0=none, 1=perturbative, 2=exact" << sampler_int << std::endl;
        exit(-1);
      }
    }
    switch (sampler_int) {
      case 0: sampler = SamplerHMC; break;
      case 1: sampler = SamplerCluster; break;
      case 2: sampler = SamplerExact; break;
      default: {
        std::cout << " ERROR: Unknown sampler: " << sampler_int << std::endl;
        std::cout << "        allowed values are 0=HMC, 1=cluster, 2=exact" << sampler_int << std::endl;
        exit(-1);
      }
    }
  } catch (std::ifstream::failure e) {      
    std::cerr << "ERROR opening or reading file \'" << filename << "\': " << e.what() << std::endl;
  }
}

/** Print out all parameters */
void Parameters::show() {
  std::string renorm_label;
  switch (renormalisation) {
    case (RenormalisationNone): renorm_label = "none"; break;
    case (RenormalisationPerturbative): renorm_label = "perturbative"; break;
    case (RenormalisationExact): renorm_label = "exact"; break;
  }
  std::string sampler_label;
  switch (sampler) {
    case (SamplerHMC): sampler_label = "HMC"; break;
    case (SamplerCluster): sampler_label = "cluster"; break;
    case (SamplerExact): sampler_label = "exact"; break;
  }
  std::cout << " M_lat            = " << M_lat << std::endl;
  std::cout << " T_final          = " << T_final << std::endl;
  std::cout << " m0               = " << m0 << std::endl;
  std::cout << " mu2              = " << mu2 << std::endl;
  std::cout << " lambda           = " << lambda << std::endl;
  std::cout << " sigma            = " << sigma << std::endl;
  std::cout << " renormalisation  = " << renorm_label << std::endl;
  std::cout << " n_burnin         = " << n_burnin << std::endl;
  std::cout << " n_samples        = " << n_samples << std::endl;
  std::cout << " sampler          = " << sampler_label << std::endl;
  std::cout << " T_hmc            = " << T_hmc << std::endl;
  std::cout << " dt_hmc           = " << dt_hmc << std::endl;
  std::cout << " n_burnin_hmc     = " << n_burnin_hmc << std::endl;
}
