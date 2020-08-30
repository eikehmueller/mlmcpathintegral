#ifndef QMPARAMETERS_HH
#define QMPARAMETERS_HH QMPARAMETERS_HH

#include <string>
#include "auxilliary/parameters.hh"

/** @file qmparameters.hh
 * 
 * @brief Parameter header file for quantum mechanical actions
 *
 */

/** @brief Enum for different actions
 *  - 0: Harmonic Oscillator
 *  - 1: Quartic Oscillator
 *  - 2: Quantum mechanical rotor
*/
enum ActionType {
  ActionHarmonicOscillator = 0,
  ActionQuarticOscillator = 1,
  ActionRotor = 2
};

/** @class QMParameters
 *
 * @brief Class for storing general quantum mechanics parameter
 * 
 */
class QMParameters : public Parameters {
public:
  /** @brief Construct a new instance */
  QMParameters() :
    Parameters("quantummechanics"),
    action_(ActionHarmonicOscillator) {
    addKey("action",String);
  }

  /** @brief Read parameters from file
   *
   * @param[in] filename Name of file to read
   */
  int readFile(const std::string filename) {

    int readSuccess = Parameters::readFile(filename);
    if (!readSuccess) {
      std::string action_str = getContents("action")->getString();
      if (action_str=="harmonicoscillator") {
        action_ = ActionHarmonicOscillator;
      } else if (action_str=="quarticoscillator") {
        action_ = ActionQuarticOscillator;
      } else if (action_str=="rotor") {
        action_ = ActionRotor;
      } else {
      }
    }
    return readSuccess;
  }

  /** @brief Return the action type */
  ActionType action() const { return action_; }

private:
  /** @brief Type of action */
  ActionType action_;
};

#endif // QMPARAMETERS_HH
