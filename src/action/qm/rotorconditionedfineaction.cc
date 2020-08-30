#include "rotorconditionedfineaction.hh"
/** @file rotorconditionedfineaction.cc
 * @brief Implementation of rotorconditionedfineaction.hh
 */


/* Fill in fine points */
void RotorConditionedFineAction::fill_fine_points(std::shared_ptr<Path> x_path) const {
  unsigned int M_lat = x_path->M_lat;
  // interior points
  for (unsigned int j=0; j<M_lat/2-1; ++j) {
    double x_m = x_path->data[2*j];
    double x_p = x_path->data[2*(j+1)];
    double x0 = action->getWminimum(x_m,x_p);
    double sigma = 2.*action->getWcurvature(x_m,x_p);
    x_path->data[2*j+1] = mod_2pi(x0 + exp_sin2_dist.draw(engine,sigma));
  }
  // Final point which requires wrap around
  double x_m = x_path->data[M_lat-2];
  double x_p = x_path->data[0];
  double x0 = action->getWminimum(x_m,x_p);
  double sigma = 2.*action->getWcurvature(x_m,x_p);
  x_path->data[M_lat-1] = mod_2pi(x0 + exp_sin2_dist.draw(engine,sigma));
}

/* Evaluate conditioned action at fine points */
double RotorConditionedFineAction::evaluate(const std::shared_ptr<Path> x_path) const {
  unsigned int M_lat = action->getM_lat();
  double x_m = x_path->data[M_lat-2];
  double x_p = x_path->data[0];
  double dx = x_path->data[M_lat-1] - action->getWminimum(x_m,x_p);
  double sigma = 2.0*action->getWcurvature(x_m,x_p);
  double S = -log(exp_sin2_dist.evaluate(dx,sigma));
  for (unsigned int j=0; j<M_lat/2-1; ++j) {
    x_m = x_path->data[2*j];
    x_p = x_path->data[2*j+2];
    double dx = x_path->data[2*j+1] - action->getWminimum(x_m,x_p);
    double sigma = 2.0*action->getWcurvature(x_m,x_p);
    S += -log(exp_sin2_dist.evaluate(dx,sigma));
  }
  return S;
}

