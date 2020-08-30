#ifndef PARAMETERS_HH
#define PARAMETERS_HH PARAMETERS_HH
#include <iostream>
#include <fstream>
#include <string>
#include <regex.h>
#include <fstream>
#include <map>
#include <vector>
#include "mpi/mpi_wrapper.hh"

/** @file parameters.hh
 * 
 * @brief Header file for parameter class
 *
 * @details Provides support for reading parameters from file and
 * storing them in appropriate classes.
 */

/** @brief Enum for different methods
 *   - 0: Single-level
 *   - 1: Two-level
 *   - 2: Multi-level
 */
enum MethodType {
  MethodSingleLevel = 0,
  MethodTwoLevel = 1,
  MethodMultiLevel = 2
};

/** @brief Enum for different actions
 *  - 0: Harmonic Oscillator
 *  - 1: Quartic Oscillator
 *  - 2: Double Well
 *  - 3: Quantum mechanical rotor
*/
enum ActionType {
  ActionHarmonicOscillator = 0,
  ActionQuarticOscillator = 1,
  ActionDoubleWell = 2,
  ActionRotor = 3
};

/** @brief Enum for sampler:
 *  - 0: HMC
 *  - 1: cluster algorithm (only some actions)
 *  - 2: exact (only some actions)
 *  - 3: hierarchical
 *  -4: multilevel
*/
enum SamplerType {
  SamplerHMC = 0,
  SamplerCluster = 1,
  SamplerExact = 2,
  SamplerHierarchical = 3,
  SamplerMultilevel = 4
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

/** @class Parameter 
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
   * @param[in] label__ Label to be used for parameter
   * @param[in] valueString String containing value
   */
  Parameter(const std::string label__,
            const std::string valueString) :
    label_(label__), valueString_(valueString) {}

  /** @brief Return value as double (if possible), return error by default */
  virtual const double getDouble() const {
    mpi_parallel::cerr << " Invalid data type (double). " << std::endl;
    return 0.0;
  };

  /** @brief Return value as integer (if possible), return error by default */
  virtual const int getInt() const {
    mpi_parallel::cerr << " Invalid data type. " << std::endl;
    return 0;
  };

  /** @brief Return value as integer (if possible), return error by default */
  virtual const std::string getString() const {
    mpi_parallel::cerr << " Invalid data type. " << std::endl;
    return std::string("");
  };

  /** @brief Return value as bool (if possible), return error by default */
  virtual const bool getBool() const {
    mpi_parallel::cerr << " Invalid data type. " << std::endl;
    return false;
  };

  /** @brief Return label */
  std::string label() const { return label_; }

  /** @brief Set new value
   *
   * @param[in] valueString String containing new value
   */
  
  // set data
  virtual void setValue(const std::string valueString) {
  };

  /** @brief Get value string */
  std::string getValueString() {
    return valueString_;
  }
  
protected:
  /** @brief string with label */
  const std::string label_;
protected:
  /** @brief string with value */
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
   * @param[in] label_ Label of parameter 
   * @param[in] valueString String containing value
   * @param[in] num_constraint_ Constraint on numerical value
   */  
  DoubleParameter(const std::string label_,
                  const std::string valueString,
                  const NumConstraintFlag num_constraint_=AnyValue) :
    Parameter(label_,valueString),
    num_constraint(num_constraint_) {
    value=atof(valueString.c_str());
  }
  /** @brief Return value as double */
  virtual const double getDouble() const { return value; }
  /** @brief Set new value
   *
   * @param[in] valueString String with value
   */
  virtual void setValue(const std::string valueString) {
    valueString_=valueString;
    value=atof(valueString.c_str());
    check_constraint();
  }

private:
  /** @brief Check constraints of current value */
  void check_constraint() {
    if ( (num_constraint == Positive) and (not (value > 0)) ) {
      mpi_parallel::cerr << "ERROR: expected positive number for parameter " << label_ << std::endl;
      mpi_exit(EXIT_FAILURE);
    }
    if ( (num_constraint == NonNegative) and (not (value >= 0)) ) {
      mpi_parallel::cerr << "ERROR: expected non-negative number for parameter " << label_ << std::endl;
      mpi_exit(EXIT_FAILURE);
    }
    if ( (num_constraint == Negative) and (not (value < 0)) ) {
      mpi_parallel::cerr << "ERROR: expected negative number for parameter " << label_ << std::endl;
      mpi_exit(EXIT_FAILURE);
    }
    if ( (num_constraint == NonPositive) and (not (value <= 0)) ) {
      mpi_parallel::cerr << "ERROR: expected non-positive number for parameter " << label_ << std::endl;
      mpi_exit(EXIT_FAILURE);
    }
  }
  /** @brief Numerical constraint flag */
  NumConstraintFlag num_constraint;
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
   * @param[in] label_ Label of parameter 
   * @param[in] valueString String containing value
   * @param[in] num_constraint Numerical constraint flag
   */  
  IntParameter(const std::string label_,
               const std::string valueString,
               const NumConstraintFlag num_constraint_=AnyValue) :
    Parameter(label_,valueString), num_constraint(num_constraint_) {
    value=atoi(valueString.c_str());
  }
  /** @brief Return value as integer */
  virtual const int getInt() const { return value; }
  /** @brief Set new value
   *
   * @param[in] valueString String with value
   */
  virtual void setValue(const std::string valueString) {
    valueString_=valueString;
    value=atoi(valueString.c_str());
    check_constraint();
  }
private:
  /** @brief Check constraints of current value */
  void check_constraint() {
    if ( (num_constraint == Positive) and (not (value > 0)) ) {
      mpi_parallel::cerr << "ERROR: expected positive number for parameter " << label_ << std::endl;
      mpi_exit(EXIT_FAILURE);
    }
    if ( (num_constraint == NonNegative) and (not (value >= 0)) ) {
      mpi_parallel::cerr << "ERROR: expected non-negative number for parameter " << label_ << std::endl;
      mpi_exit(EXIT_FAILURE);
    }
    if ( (num_constraint == Negative) and (not (value < 0)) ) {
      mpi_parallel::cerr << "ERROR: expected negative number for parameter " << label_ << std::endl;
      mpi_exit(EXIT_FAILURE);
    }
    if ( (num_constraint == NonPositive) and (not (value <= 0)) ) {
      mpi_parallel::cerr << "ERROR: expected non-positive number for parameter " << label_ << std::endl;
      mpi_exit(EXIT_FAILURE);
    }
  }
  /** @brief Numerical constraint flag */
  NumConstraintFlag num_constraint;
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
   * @param[in] label_ Label of parameter 
   * @param[in] valueString String containing value
   */  
  StringParameter(const std::string label_,
                  const std::string valueString) :
    Parameter(label_,valueString) {
    value = valueString.substr(1, valueString.size()-2);
  }
  /** @brief Return value as string */
  virtual const std::string getString() const { return value; }
  /** @brief Set new value
   *
   * @param[in] valueString String with value
   */
  virtual void setValue(const std::string valueString) {
    valueString_=valueString;
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
   * @param[in] label_ Label of parameter 
   * @param[in] valueString String containing value
   */  
  BoolParameter(const std::string label,
                const std::string valueString) :
    Parameter(label,valueString) {
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
    valueString_=valueString;
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
   * @param[in] num_constraint Numerical constraint flag
   */
  void addKey(const std::string key,
              const Datatype datatype,
              const NumConstraintFlag num_constraint=AnyValue);
  /** @brief hash for key - value pairs */
  std::map<std::string,std::vector<Parameter*> > contents;
private:
  /** @brief hash for key - datatype pairs */
  std::map<std::string,Datatype> keywords;
  /** @brief hash for numerical constraint flags */
  std::map<std::string,NumConstraintFlag> num_constraints;
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

/** @class GeneralParameters
 *
 * @brief Class for storing general parameter
 * 
 */
class GeneralParameters : public Parameters {
public:
  /** @brief Construct a new instance */
  GeneralParameters() :
    Parameters("general"),
    method_(MethodSingleLevel),
    action_(ActionHarmonicOscillator) {
    addKey("method",String);
    addKey("action",String);
  }

  /** @brief Read parameters from file
   *
   * @param[in] filename Name of file to read
   */
  int readFile(const std::string filename) {

    int readSuccess = Parameters::readFile(filename);
    if (!readSuccess) {
      std::string method_str = getContents("method")->getString();
      if (method_str=="singlelevel") {
        method_ = MethodSingleLevel;
      } else if (method_str=="twolevel") {
        method_ = MethodTwoLevel;
      } else if (method_str=="multilevel") {
        method_ = MethodMultiLevel;
      } else {
      }
      std::string action_str = getContents("action")->getString();
      if (action_str=="harmonicoscillator") {
        action_ = ActionHarmonicOscillator;
      } else if (action_str=="quarticoscillator") {
        action_ = ActionQuarticOscillator;
      } else if (action_str=="doublewell") {
        action_ = ActionDoubleWell;
      } else if (action_str=="rotor") {
        action_ = ActionRotor;
      } else {
      }
    }
    return readSuccess;
  }

  /** @brief Method to use (single-, two-, or multi-level  */
  MethodType method() const { return method_; }
  /** @brief Return the action type */
  ActionType action() const { return action_; }

private:
  /** @brief Method to use */
  MethodType method_;
  /** @brief Type of action */
  ActionType action_;
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
    addKey("M_lat",Integer,Positive);
    addKey("T_final",Double,Positive);
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

/** @brief Write parameters to stream */
std::ostream& operator<<(std::ostream& output, const Parameters& p);

#endif // PARAMETERS_HH
