#include "timer.hh"
/** @file timer.cc
 * @brief Implementation of timer.hh
 */

/* Start timer */
void Timer::start() {
  mpi_barrier();
  tp_start = Clock::now();
}

/* Stop timer */
void Timer::stop() {
  mpi_barrier();
  tp_finish = Clock::now();
  std::chrono::duration<double> t_duration = tp_finish - tp_start;
  t_elapsed = double(t_duration.count());
}

/* Reset timer */
void Timer::reset() {
  stop();
  t_elapsed = 0.0;
}

/* Write to stream */
std::ostream &operator<<(std::ostream &out, const Timer &t) {
  out.precision(4);
  out << "[timer " << t.Label() << "] : " << std::scientific << t.elapsed()
      << " s";
  return out;
}

/* return current time as string */
std::string current_time() {
  auto time_now = std::chrono::system_clock::now();
  std::time_t t = std::chrono::system_clock::to_time_t(time_now);
  return std::ctime(&t);
}
