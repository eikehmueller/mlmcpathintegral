#ifndef TIMER_HH
#define TIMER_HH TIMER_HH

#include "mpi/mpi_wrapper.hh"
#include <chrono>
#include <ctime>
#include <iostream>
#include <string>
#include <time.h>

/** @file timer.hh
 * @brief Header file for timer class
 */

/** @class Timer
 *
 * @brief Timer class
 *
 * This class allows simple timer measurements using the chrono library
 */

/* ************************************* *
 * Class for timer
 * ************************************* */
class Timer {
public:
  typedef std::chrono::high_resolution_clock Clock;
  typedef std::chrono::time_point<Clock> Timepoint;
  /** @brief Construct rew instance
   *
   * @param[in] label Label for identifying timer in output
   */
  Timer(const std::string label_ = "") : label(label_), t_elapsed(0.0) {}
  /** @brief Return elapsed timer in seconds */
  double elapsed() const { return t_elapsed; }
  /** @brief Reset timer */
  void reset();
  /** @brief Start timer again */
  void start();
  /** @brief Stop timer */
  void stop();
  /** @brief Return label string */
  std::string Label() const { return label; }

private:
  /** @brief Unique label string */
  std::string label;
  /** @brief total elapsed time */
  double t_elapsed;
  /** @brief Internal start timepoint */
  Timepoint tp_start;
  /** @brief Internal end timepoint */
  Timepoint tp_finish;
};

/** @brief Current time
 * Return current time as a string
 */
std::string current_time();

/** @brief Write timer to stream
 *
 * @param[in] t Instance of timer class to print out
 * @param[inout] out Output stream to write to
 */
std::ostream &operator<<(std::ostream &out, const Timer &t);

#endif // TIMER_HH
