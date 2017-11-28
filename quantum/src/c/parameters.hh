#ifndef PARAMETERS_HH
#define PARAMETERS_HH PARAMETERS_HH
#include <iostream>
#include <fstream>
#include <string>

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
   * M_lat     VALUE
   * T_final   VALUE
   * m0        VALUE
   * mu2       VALUE
   * n_burnin  VALUE
   * n_samples VALUE
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
  std::string filename; /** Name of file to read */
  unsigned int M_lat; /** Number of time slices */
  double T_final; /** Final time */
  double m0; /** Particle mass */
  double mu2; /** parameter \f$\mu^2\f$ in harmonic oscillator potential */
  unsigned int n_burnin; /** Number of burnin-samples */
  unsigned int n_samples; /** Number of samples */
};
#endif // PARAMETERS_HH
