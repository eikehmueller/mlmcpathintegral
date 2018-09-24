#include "timer.hh"
/** @file timer.cc
 * @brief Implementation of timer.hh
 */


/* Start timer */
void Timer::start() {
  tp_start = Clock::now();
}

/* Stop timer */
void Timer::stop() {
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
std::ostream& operator<<(std::ostream& out, const Timer& t) {
  out.precision(4);
  out << "[timer " << t.Label() << "] : " << std::scientific << t.elapsed() << " s";
  return out;
}
