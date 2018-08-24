#ifndef PARAMETERS_HH
#define PARAMETERS_HH PARAMETERS_HH
#include <iostream>
#include <fstream>
#include <string>

/** @file parameters.hh
 * @brief Header file for parameter class
 */

/** @brief Enum for renormalisations:
 *  - 0: No renormalisation
 *  - 1: Perturbative renormalisation
 *  - 2: Exact renormalisation
*/
enum RenormalisationType {
  RenormalisationNone = 0,
  RenormalisationPerturbative = 1,
  RenormalisationExact = 2
};

/** @brief Enum for sampler:
 *  - 0: HMC
 *  - 1: cluster algorithm (only some actions)
 *  - 2: exact (only some actions)
*/
enum SampleType {
  SamplerHMC = 0,
  SamplerCluster = 1,
  SamplerExact = 2
};

/** @class Parameters
 *
 * @brief Storage for parameters, which can be read from a file 
 */
class Parameters {
public:
  /** @brief Create new instance
   *
   * Read parameters from specified file. The content of the file has to be:
   *
   @verbatim
     M_lat           VALUE
     T_final         VALUE
     m0              VALUE
     mu2             VALUE
     lambda          VALUE
     sigma           VALUE
     renormalisation VALUE
     n_burnin        VALUE
     n_samples       VALUE
     sampler         VALUE
     T_hmc           VALUE
     dt_hmc          VALUE
     n_burnin_hmc    VALUE
   @endverbatim
   *
   * @param[in] filename_ Name of file to read
   */
  Parameters(const std::string filename_) : filename(filename_) {
    readfile();
  }

  /** @brief Read parameters from file */
  void readfile();

  /** @brief Print all parameters to screen */
  void show();

public:
  /** @brief Name of file to read */
  std::string filename;
  /** @brief Number of time slices */
  unsigned int M_lat;
  /** @brief Final time */
  double T_final;
  /** @brief Particle mass */
  double m0;
  /** @brief parameter \f$\mu^2\f$ in harmonic oscillator potential */
  double mu2;
  /** @brief parameter \f$\lambda\f$ in quartic oscillator potential */
  double lambda;
  /** @brief width parameter \f$\sigma\f$ in double well potential */
  double sigma;
  /** @brief Renormalisation type (harmonic oscillator only) */
  RenormalisationType renormalisation;
  /** @brief Number of burn-in samples */
  unsigned int n_burnin;
  /** @brief Number of samples */
  unsigned int n_samples;
  /** @brief Use Hybrid Monte Carlo (HMC) sampling? */
  int sampler;
  /** @brief Length \f$T_{hmc}\f$ of HMC trajectories */
  double T_hmc;
  /** @brief Time step size \f$\Delta t_{hmc}\f$ of HMC trajectories */
  double dt_hmc;
  /** @brief HMC burnin samples */
  unsigned int n_burnin_hmc;
};
#endif // PARAMETERS_HH
