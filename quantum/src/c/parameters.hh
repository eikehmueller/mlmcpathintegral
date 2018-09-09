#ifndef PARAMETERS_HH
#define PARAMETERS_HH PARAMETERS_HH
#include <iostream>
#include <fstream>
#include <string>
#include <regex.h>
#include<fstream>
#include<map>
#include<vector>

/** @file parameters.hh
 * 
 * @brief Header file for parameter class
 *
 * @details Provides support for reading parameters from file and
 * storing them in appropriate classes.
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
enum SamplerType {
  SamplerHMC = 0,
  SamplerCluster = 1,
  SamplerExact = 2
};

/** @brief Flags with constraint on numbers
 * @details Allows restricting values of parameters
 */
enum NumConstraintFlag {
  AnyValue = 0,    // x can be any value
  Positive = 1,    // x > 0
  NonNegative = 2, // x >= 0
  Negative = 3,    // x < 0
  NonPositive = 4  // x <=0
};

/** @class Paramater 
 *
 *  @brief Base class for storing a parameter
 * 
 *  @details This class is later overloaded for parameters of different
 *           types.
 */
class Parameter {
public:
  /** @brief Construct a new instance
   *
   * @param[in] valueString String containing value
   */
  Parameter(const std::string valueString) :
    valueString_(valueString) {}

  /** @brief Return value as double (if possible), return error by default */
  virtual const double getDouble() const {
    std::cerr << " Invalid data type (double). " << std::endl;
    return 0.0;
  };

  /** @brief Return value as integer (if possible), return error by default */
  virtual const int getInt() const {
    std::cerr << " Invalid data type. " << std::endl;
    return 0;
  };

  /** @brief Return value as integer (if possible), return error by default */
  virtual const std::string getString() const {
    std::cerr << " Invalid data type. " << std::endl;
    return std::string("");
  };

  /** @brief Return value as bool (if possible), return error by default */
  virtual const bool getBool() const {
    std::cerr << " Invalid data type. " << std::endl;
    return false;
  };

  /** @brief Set new value
   *
   * @param[in] valueString String containing new value
   */
  // set data
  virtual void setValue(const std::string valueString) {
  };
private:
  std::string valueString_;
};

/** @class DoubleParameter 
 *
 * @brief Class for representing a double parameter
 */
class DoubleParameter : public Parameter {
public:
  /** @brief Construct a new instance
   *
   * @param[in] valueString String containing value
   */  
  DoubleParameter(const std::string valueString,
                  const NumConstraintFlag numconstraint=AnyValue) :
    Parameter(valueString) {
    value=atof(valueString.c_str());
    if ( (numconstraint == Positive) and (not (value > 0)) ) {
      std::cerr << "ERROR: expected positive number" << std::endl;
    }
    if ( (numconstraint == NonNegative) and (not (value >= 0)) ) {
      std::cerr << "ERROR: expected non-negative number" << std::endl;
    }
    if ( (numconstraint == Negative) and (not (value < 0)) ) {
      std::cerr << "ERROR: expected negative number" << std::endl;
    }
    if ( (numconstraint == NonPositive) and (not (value <= 0)) ) {
      std::cerr << "ERROR: expected non-positive number" << std::endl;
    }
  }
  /** @brief Return value as double */
  virtual const double getDouble() const { return value; }
  /** @brief Set new value
   *
   * @param[in] valueString String with value
   */
  virtual void setValue(const std::string valueString) {
    value=atof(valueString.c_str());
  }
private:
  /** @brief Value of parameter */
  double value;
};

/** @class IntParameter 
 *
 * @brief Class for representing an integer parameter
 */
class IntParameter : public Parameter {
public:
  /** @brief Construct a new instance
   *
   * @param[in] valueString String containing value
   */  
  IntParameter(const std::string valueString,
               const NumConstraintFlag numconstraint=AnyValue) :
    Parameter(valueString) {
    value=atoi(valueString.c_str());
    if ( (numconstraint == Positive) and (not (value > 0)) ) {
      std::cerr << "ERROR: expected positive number" << std::endl;
    }
    if ( (numconstraint == NonNegative) and (not (value >= 0)) ) {
      std::cerr << "ERROR: expected non-negative number" << std::endl;
    }
    if ( (numconstraint == Negative) and (not (value < 0)) ) {
      std::cerr << "ERROR: expected negative number" << std::endl;
    }
    if ( (numconstraint == NonPositive) and (not (value <= 0)) ) {
      std::cerr << "ERROR: expected non-positive number" << std::endl;
    }
  }
  /** @brief Return value as integer */
  virtual const int getInt() const { return value; }
  /** @brief Set new value
   *
   * @param[in] valueString String with value
   */
  virtual void setValue(const std::string valueString) {
    value=atoi(valueString.c_str());
  }
private:
  /** @brief Value of parameter */
  int value;
};

/** @class StringParameter 
 *
 * @brief Class for representing a string valued parameter
 */
class StringParameter : public Parameter {
public:
  /** @brief Construct a new instance
   *
   * @param[in] valueString String containing value
   */  
  StringParameter(const std::string valueString) :
    Parameter(valueString) {
    value = valueString.substr(1, valueString.size()-2);
  }
  /** @brief Return value as string */
  virtual const std::string getString() const { return value; }
  /** @brief Set new value
   *
   * @param[in] valueString String with value
   */
  virtual void setValue(const std::string valueString) {
    value = valueString.substr(1, valueString.size()-2);
  }
private:
  /** @brief Value of parameter */
  std::string value;
};

/** @class DoubleParameter 
 *
 * @brief Class for representing a double parameter
 */
class BoolParameter : public Parameter {
public:
  /** @brief Construct a new instance
   *
   * @param[in] valueString String containing value
   */  
  BoolParameter(const std::string valueString) :
    Parameter(valueString) {
    value = ( ( valueString == "true" ) ||
              ( valueString == "True" ) ||
              ( valueString == "TRUE" ) ) ;

  }
  /** @brief Return value as double */
  virtual const bool getBool() const { return value; }
  /** @brief Set new value
   *
   * @param[in] valueString String with value
   */
  virtual void setValue(const std::string valueString) {
    value = ( ( valueString == "true" ) ||
              ( valueString == "True" ) ||
              ( valueString == "TRUE" ) ) ;
  }
private:
  /** @brief Value of parameter */
  bool value;
};

/** @class Parameters
 *
 * @brief Base class for parameters which can be read from a file.
 *
 * @details The file has to have the following structure:
 *
 * @verbatim
 *
 * sectionname:
 *   key = value
 *   key = value1, value2
 *   ...
 * sectionname:
 *   key = value
 *   key = value1, value2, value3
 *   ...
 * @endverbatim
 *
 * The section to be parsed by this particular class is given by the
 * sectionname which is passed to the constructor. Blank lines and lines
 * starting with '#' are ignored. Comments starting with '#' at the end of
 * a line are also ignored.
 *
 * ************************************************************** */
class Parameters {
public:
  enum Datatype {Integer, Double, String, Bool};
  /** @brief allow access for output redirection */
  friend std::ostream& operator<<(std::ostream& output, const Parameters& p);
  
  /** @brief create a new instance
   * 
   * @param[in] section Name of the section to be parsed
   */
  Parameters(const std::string section);

  /** @brief Read parameters from file
   * 
   * This reads all data from the file and extracts the data into a form
   * which can later be extracted with the getContents() methods.
   * Internally the data is stored in a map object.
   * 
   * @param[in] filename Name of file to be read
   */
  int readFile(const std::string filename);

  /** @brief Extract scalar data associated with a particular key 
   * 
   * This returns the data as a Parameter object.
   *
   * @param[in] key Name of key to be read
   */
  const Parameter* getContents(const std::string key) {
    return contents[key][0];
  }

  /** @brief Extract vector data associated with a particular key 
   * 
   * This returns the data as a Parameter object.
   *
   * @param[in] key Name of key to be read
   */
  const Parameter* getContents(const std::string key, const int i) {
    return contents[key][i];
  }

private:
  /** @brief Auxilliary method for splitting a string and removing whitespace
   *
   * Split a string by delimiter ',' and remove all whitespace.
   *
   * @param[in] s String to be processes
   */
  std::vector<std::string> split_string(std::string s);

protected:
  /** @brief Register a particular key-value pair for specific data type 
   * 
   * @param[in] key Name of key to be used
   * @param[in] datatype Datatype of value
   */
  void addKey(const std::string key, const Datatype datatype);
  /** @brief hash for key - value pairs */
  std::map<std::string,std::vector<Parameter*> > contents;
private:
  /** @brief hash for key - datatype pairs */
  std::map<std::string,Datatype> keywords;
  /** @brief Name of section associated with this class */
  std::string section_;
  /** @brief Regular expression for a comment */
  regex_t regexComment_;
  /** @brief Regular expression for a keyword */
  regex_t regexKeyword_;
  /** @brief Regular expression for a key-value pair */
  regex_t regexKeyValueAny_;
  /** @brief Regular expression for an integer value */
  regex_t regexKeyValueInt_;
  /** @brief Regular expression for a double value */
  regex_t regexKeyValueDouble_;
  /** @brief Regular expression for a bool value */
  regex_t regexKeyValueBool_;
  /** @brief Regular expression for a string value */
  regex_t regexKeyValueString_;
  /** @brief Verbosity flag */
  static const bool verbose=0;
};

/** @class GeneralParameter
 *
 * @brief Class for storing general parameter
 * 
 */
class GeneralParameters : public Parameters {
public:
  /** @brief Construct a new instance */
  GeneralParameters() :
    Parameters("general"),
    do_singlelevelmc_(true),
    do_twolevelmc_(true) {
    addKey("do_singlelevelmc",Bool);
    addKey("do_twolevelmc",Bool);
  }

  /** @brief Read parameters from file
   *
   * @param[in] filename Name of file to read
   */
  int readFile(const std::string filename) {

    int readSuccess = Parameters::readFile(filename);
    if (!readSuccess) {
      do_singlelevelmc_ = getContents("do_singlelevelmc")->getBool();
      do_twolevelmc_ = getContents("do_twolevelmc")->getBool();
    }
    return readSuccess;
  }

  /** @brief Run single level MC algorithm? */
  bool do_singlelevelmc() const { return do_singlelevelmc_; }
  /** @brief Run two level MC algorithm? */
  bool do_twolevelmc() const { return do_twolevelmc_; }

private:
  /** @brief Run single level MC algorithm? */
  bool do_singlelevelmc_;
  /** @brief Run two level MC algorithm? */
  bool do_twolevelmc_;
};


/** @class LatticeParameters
 *
 * @brief Class for storing parameters of lattice
 * 
 * This stores the number \f$M_{lat}\f$ of lattice sites and the
 * final time \f$T_{final}\f$
 */
class LatticeParameters : public Parameters {
public:
  /** @brief Construct a new instance */
  LatticeParameters() :
    Parameters("lattice"),
    M_lat_(8),
    T_final_(1.0) {
    addKey("M_lat",Integer);
    addKey("T_final",Double);
  }

  /** @brief Read parameters from file
   *
   * @param[in] filename Name of file to read
   */
  int readFile(const std::string filename) {

    int readSuccess = Parameters::readFile(filename);
    if (!readSuccess) {
      M_lat_ = getContents("M_lat")->getInt();
      T_final_ = getContents("T_final")->getDouble();
    }
    return readSuccess;
  }

  /** @brief Return number of lattice sites */
  unsigned int M_lat() const { return M_lat_; }
  /** @brief Return final time */
  double T_final() const { return T_final_; }

private:
  /** @brief Number of lattice sites */
  unsigned int M_lat_;
  /** @brief Final time */
  double T_final_;
};

/** @class HarmonicOscillatorParameters
 *
 * @brief Class for storing parameters of harmonic oscillator action
 *
 * This stores the mass \f$m_0\f$ and curvature \f$\mu_2\f$ of the 
 * harmonic oscillator action with potential \f$V(x)=\frac{m_0}{2}\mu^2x^2\f$
 */
class HarmonicOscillatorParameters : public Parameters {
public:
  /** @brief Construct a new instance */
  HarmonicOscillatorParameters() :
    Parameters("harmonicoscillator"),
    m0_(1.0),
    mu2_(1.0),
    renormalisation_(RenormalisationNone) {
    addKey("m0",Double);
    addKey("mu2",Double);
    addKey("renormalisation",String);
  }

  /** @brief Read parameters from file
   *
   * @param[in] filename Name of file to read
   */
  int readFile(const std::string filename) {

    int readSuccess = Parameters::readFile(filename);
    if (!readSuccess) {
      m0_ = getContents("m0")->getDouble();
      mu2_ = getContents("mu2")->getDouble();
      std::string renormalisation_str = getContents("renormalisation")->getString();
      if (renormalisation_str == "none") {
        renormalisation_ = RenormalisationNone;
      } else if (renormalisation_str == "perturbative") {
        renormalisation_ = RenormalisationPerturbative;
      } else if (renormalisation_str == "exact") {
        renormalisation_ = RenormalisationExact;
      }
    }
    return readSuccess;
  }

  /** @brief Return unrenormalised mass \f$m_0\f$ */
  double m0() const { return m0_; }
  /** @brief Return parameter \f$\mu^2\f$ */
  double mu2() const { return mu2_; }
  /** @brief Return renormalisation */
  RenormalisationType renormalisation() const { return renormalisation_; }

private:
  /** @brief Unrenormalised mass \f$m_0\f$ */
  double m0_;
  /** @brief Parameter \f$\mu^2\f$ */
  double mu2_;
  /** @brief Renormalisation */
  RenormalisationType renormalisation_;
};

/** @class DoubleWellParameters
 *
 * @brief Class for storing parameters of double well action
 *
 * This stores the mass \f$m_0\f$ and parameters \f$\mu_2\f$, \f$\lambda\f$
 * and \f$\sigma\f$ of the double well action with potential
 * \f$V(x)=\frac{m_0}{2}\mu^2x^2+\lambda\exp\left(-\frac{x^2}{2\sigma^2}\right)\f$
 */
class DoubleWellParameters : public Parameters {
public:
  /** @brief Construct a new instance */
  DoubleWellParameters() :
    Parameters("doublewell"),
    m0_(1.0),
    mu2_(1.0),
    lambda_(1.0),
    sigma_(1.0) {
    addKey("m0",Double);
    addKey("mu2",Double);
    addKey("lambda",Double);
    addKey("sigma",Double);
  }

  /** @brief Read parameters from file
   *
   * @param[in] filename Name of file to read
   */
  int readFile(const std::string filename) {

    int readSuccess = Parameters::readFile(filename);
    if (!readSuccess) {
      m0_ = getContents("m0")->getDouble();
      mu2_ = getContents("mu2")->getDouble();
      lambda_ = getContents("lambda")->getDouble();
      sigma_ = getContents("sigma")->getDouble();
    }
    return readSuccess;
  }

  /** @brief Return unrenormalised mass \f$m_0\f$ */
  double m0() const { return m0_; }
  /** @brief Return parameter \f$\mu^2\f$ */
  double mu2() const { return mu2_; }
  /** @brief Return parameter \f$\lambda\f$ */
  double lambda() const { return lambda_; }
  /** @brief Return parameter \f$\sigma\f$ */
  double sigma() const { return sigma_; }

private:
  /** @brief Unrenormalised mass \f$m_0\f$ */
  double m0_;
  /** @brief Parameter \f$\mu^2\f$ */
  double mu2_;
  /** @brief Parameter \f$\lambda\f$ */
  double lambda_;
  /** @brief Parameter \f$\sigma\f$ */
  double sigma_;
};

/** @class QuarticOscillatorParameters
 *
 * @brief Class for storing parameters of quartic oscillator action
 *
 * This stores the mass \f$m_0\f$ and parameters \f$\mu_2\f$, \f$\lambda\f$
 * of the double well action with potential
 * \f$V(x)=\frac{m_0}{2}\mu^2x^2+\frac{\lambda}{4}x^4\f$.
 */
class QuarticOscillatorParameters : public Parameters {
public:
  /** @brief Construct a new instance */
  QuarticOscillatorParameters() :
    Parameters("quarticoscillator"),
    m0_(1.0),
    mu2_(1.0),
    lambda_(1.0) {
    addKey("m0",Double);
    addKey("mu2",Double);
    addKey("lambda",Double);
  }

  /** @brief Read parameters from file
   *
   * @param[in] filename Name of file to read
   */
  int readFile(const std::string filename) {

    int readSuccess = Parameters::readFile(filename);
    if (!readSuccess) {
      m0_ = getContents("m0")->getDouble();
      mu2_ = getContents("mu2")->getDouble();
      lambda_ = getContents("lambda")->getDouble();
    }
    return readSuccess;
  }

  /** @brief Return unrenormalised mass \f$m_0\f$ */
  double m0() const { return m0_; }
  /** @brief Return parameter \f$\mu^2\f$ */
  double mu2() const { return mu2_; }
  /** @brief Return parameter \f$\lambda\f$ */
  double lambda() const { return lambda_; }

private:
  /** @brief Unrenormalised mass \f$m_0\f$ */
  double m0_;
  /** @brief Parameter \f$\mu^2\f$ */
  double mu2_;
  /** @brief Parameter \f$\lambda\f$ */
  double lambda_;
};

/** @class RotorParameters
 *
 * @brief Class for storing parameters of quantum rotor action
 *
 * This stores the mass \f$m_0\f$ of the quantum mechanical rotor.
 */
class RotorParameters : public Parameters {
public:
  /** @brief Construct a new instance */
  RotorParameters() :
    Parameters("rotor"),
    m0_(1.0) {
    addKey("m0",Double);
  }

  /** @brief Read parameters from file
   *
   * @param[in] filename Name of file to read
   */
  int readFile(const std::string filename) {

    int readSuccess = Parameters::readFile(filename);
    if (!readSuccess) {
      m0_ = getContents("m0")->getDouble();
    }
    return readSuccess;
  }

  /** @brief Return unrenormalised mass \f$m_0\f$ */
  double m0() const { return m0_; }

private:
  /** @brief Unrenormalised mass \f$m_0\f$ */
  double m0_;
};

/** @class SingleLevelMCParameters
 *
 * @brief Class for storing parameters of single level Monte Carlo integrator.
 */
class SingleLevelMCParameters : public Parameters {
public:
  /** @brief Construct a new instance */
  SingleLevelMCParameters() :
    Parameters("singlelevelmc"),
    n_burnin_(100),
    n_samples_(100),
    sampler_(SamplerHMC) {
    addKey("n_burnin",Integer);
    addKey("n_samples",Integer);
    addKey("sampler",String);
  }

  /** @brief Read parameters from file
   *
   * @param[in] filename Name of file to read
   */
  int readFile(const std::string filename) {

    int readSuccess = Parameters::readFile(filename);
    if (!readSuccess) {
      n_burnin_ = getContents("n_burnin")->getInt();
      n_samples_ = getContents("n_samples")->getInt();
      std::string sampler_str = getContents("sampler")->getString();
      if (sampler_str == "HMC") {
        sampler_ = SamplerHMC;
      } else if (sampler_str == "cluster") {
        sampler_ = SamplerCluster;
      } else if (sampler_str == "exact") {
        sampler_ = SamplerExact;
      } else  {
        std::cerr << " ERROR: Unknown sampler: " << sampler_str << std::endl;
        std::cerr << "        allowed values are \'HMC\', \'cluster\', \'exact\'" << std::endl;
        exit(-1);
      }
    }
    return readSuccess;
  }

  /** @brief Return number of burnin samples */
  unsigned int n_burnin() const { return n_burnin_; }
  /** @brief Return number of samples */
  unsigned int n_samples() const { return n_samples_; }
  /** @brief Return sampler type */
  SamplerType sampler() const { return sampler_; }
private:
  /** @brief Number of burnin samples */
  unsigned int n_burnin_;
  /** @brief Number of samples */
  unsigned int n_samples_;
  /** @brief Sampler type */
  SamplerType sampler_;
};

/** @class TwoLevelMCParameters
 *
 * @brief Class for storing parameters of two level Monte Carlo integrator.
 */
class TwoLevelMCParameters : public Parameters {
public:
  /** @brief Construct a new instance */
  TwoLevelMCParameters() :
    Parameters("twolevelmc"),
    n_burnin_(100),
    n_samples_(100),
    coarsesampler_(SamplerHMC) {
    addKey("n_burnin",Integer);
    addKey("n_samples",Integer);
    addKey("coarsesampler",String);
  }

  /** @brief Read parameters from file
   *
   * @param[in] filename Name of file to read
   */
  int readFile(const std::string filename) {

    int readSuccess = Parameters::readFile(filename);
    if (!readSuccess) {
      n_burnin_ = getContents("n_burnin")->getInt();
      n_samples_ = getContents("n_samples")->getInt();
      std::string sampler_str = getContents("coarsesampler")->getString();
      if (sampler_str == "HMC") {
        coarsesampler_ = SamplerHMC;
      } else if (sampler_str == "cluster") {
        coarsesampler_ = SamplerCluster;
      } else if (sampler_str == "exact") {
        coarsesampler_ = SamplerExact;
      } else  {
        std::cerr << " ERROR: Unknown coarse sampler: " << sampler_str;
        std::cerr << std::endl;
        std::cerr << "        allowed values are \'HMC\', \'cluster\', \'exact\'" << std::endl;
        exit(-1);
      }
    }
    return readSuccess;
  }

  /** @brief Return number of burnin samples */
  unsigned int n_burnin() const { return n_burnin_; }
  /** @brief Return number of samples */
  unsigned int n_samples() const { return n_samples_; }
  /** @brief Return sampler type */
  SamplerType coarsesampler() const { return coarsesampler_; }
private:
  /** @brief Number of burnin samples */
  unsigned int n_burnin_;
  /** @brief Number of samples */
  unsigned int n_samples_;
  /** @brief Sampler type */
  SamplerType coarsesampler_;
};

/** @class HMCParameters
 *
 * @brief Class for storing parameters of Hybrid MC sampler
 */
class HMCParameters : public Parameters {
public:
  /** @brief Construct a new instance */
  HMCParameters() :
    Parameters("hmc"),
    T_(1.0),
    dt_(0.1),
    n_burnin_(100) {
    addKey("T",Double);
    addKey("dt",Double);
    addKey("n_burnin",Integer);
  }

  /** @brief Read parameters from file
   *
   * @param[in] filename Name of file to read
   */
  int readFile(const std::string filename) {

    int readSuccess = Parameters::readFile(filename);
    if (!readSuccess) {
      T_ = getContents("T")->getDouble();
      dt_ = getContents("dt")->getDouble();
      n_burnin_ = getContents("n_burnin")->getInt();
    }
    return readSuccess;
  }
  /** @brief Return length of integration interval \f$T\f$ */
  double T() const { return T_; }
  /** @brief Return integration timestep \f$dt\f$ */
  double dt() const { return dt_; }
  /** @brief Return number of burnin samples */
  unsigned int n_burnin() const { return n_burnin_; }
private:
  /** @brief Integration time interval \f$T\f$ */
  double T_;
  /** @brief Integration time step \f$dt\f$ */
  double dt_;
  /** @brief Number of burnin samples */
  unsigned int n_burnin_;
};

/** @class ClusterParameters
 *
 * @brief Class for storing parameters of Cluster algorithm
*/
class ClusterParameters : public Parameters {
public:
  /** @brief Construct a new instance */
  ClusterParameters() :
    Parameters("clusteralgorithm"),
    n_burnin_(100) {
    addKey("n_burnin",Integer);
  }

  /** @brief Read parameters from file
   *
   * @param[in] filename Name of file to read
   */
  int readFile(const std::string filename) {

    int readSuccess = Parameters::readFile(filename);
    if (!readSuccess) {
      n_burnin_ = getContents("n_burnin")->getInt();
    }
    return readSuccess;
  }
  /** @brief Return number of burnin samples */
  unsigned int n_burnin() const { return n_burnin_; }
private:
  /** @brief Number of burnin samples */
  unsigned int n_burnin_;
};

/** @brief Write parameters to stream */
std::ostream& operator<<(std::ostream& output, const Parameters& p);

#endif // PARAMETERS_HH
