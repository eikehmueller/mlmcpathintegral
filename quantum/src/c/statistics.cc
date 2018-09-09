#include "statistics.hh"

std::ostream& operator<<(std::ostream& os, const Statistics& stats) {
  os << " ";
  os << std::setprecision(6) << std::fixed;
  os << stats.label() << ": Avg +/- Err = " << stats.average();
  os << " +/- " << stats.error() << std::endl;
  os << " " << stats.label() << ": Var         = " << stats.variance() << std::endl;
  os << std::setprecision(3) << std::fixed;
  os << " " << stats.label() << ": tau_{int}   = " << stats.tau_int() << std::endl;
  return os;
}
