#include "rotoraction.hh"
#include <iostream>
/** @file rotoraction.cc
 * @brief Implementation of rotoraction.hh
 */

/** Evaluate action */
const double RotorAction::evaluate(const Path* x_path) const {
  double x_diff = x_path->data[0]-x_path->data[M_lat-1];
  double S = 1.-cos(x_diff);
  for (unsigned int j=1;j<M_lat;++j) {
    double x_diff = x_path->data[j]-x_path->data[j-1];
    S += 1.-cos(x_diff);
  }
  return m0/a_lat*S;
}

/** Calculate force */
void RotorAction::force(const Path* x_path,
                        Path* p_path) const {
  double tmp = m0/a_lat;
  // Left boundary
  double x_m = x_path->data[M_lat-1];
  double x = x_path->data[0];
  double x_p = x_path->data[1];
  p_path->data[0] = tmp*(sin(x-x_m)+sin(x-x_p));
  // Interior points
  for (unsigned int j=1;j<M_lat-1;++j) {
    x_m = x_path->data[j-1];
    x = x_path->data[j];
    x_p = x_path->data[j+1];  
    p_path->data[j] = tmp*(sin(x-x_m)+sin(x-x_p));
  }
  // Right boundary
  x_m = x_path->data[M_lat-2];
  x = x_path->data[M_lat-1];
  x_p = x_path->data[0];
  p_path->data[M_lat-1] = tmp*(sin(x-x_m)+sin(x-x_p));
}

/** Initialise path */
void RotorAction::initialise_path(Path* x_path) const {
  std::mt19937_64 engine;
  double pi = 4.0*atan(1.0);
  std::uniform_real_distribution<double> uniform(-pi,pi);
  std::generate(x_path->data,x_path->data+M_lat,[&uniform,&engine]() {return uniform(engine);});
  x_path->save_to_disk("path_initial.dat");
}