#include "conditionedfineaction.hh"
/** @file conditionedfineaction.cc
 * @brief Implementation of conditionedfineaction.hh
 */

/* Fill in fine points */
void GaussianConditionedFineAction::fill_fine_points(Path* x_path) const {
  unsigned int M_lat = x_path->M_lat;
  // interior points
  for (unsigned int j=0; j<M_lat/2-1; ++j) {
    double x_m = x_path->data[2*j];
    double x_p = x_path->data[2*(j+1)];
    double x0 = action.getWminimum(x_m,x_p);
    double sigma = 1./sqrt(action.getWcurvature(x_m,x_p));
    x_path->data[2*j+1] = x0 + normal_dist(engine)*sigma;
  }
  // Final point which requires wrap around
  double x_m = x_path->data[M_lat-2];
  double x_p = x_path->data[0];
  double x0 = action.getWminimum(x_m,x_p);
  double sigma = 1./sqrt(action.getWcurvature(x_m,x_p));
  x_path->data[M_lat-1] = x0 + normal_dist(engine)*sigma;
}

/* Evaluate conditioned action at fine points */
double GaussianConditionedFineAction::evaluate(const Path* x_path) const {
  unsigned int M_lat = action.getM_lat();
  double x_m = x_path->data[M_lat-2];
  double x_p = x_path->data[0];
  double dx = x_path->data[M_lat-1] - action.getWminimum(x_m,x_p);
  double curvature = action.getWcurvature(x_m,x_p);
  double S = 0.5*curvature*dx*dx + 0.5*log(curvature);
  for (unsigned int j=0; j<M_lat/2-1; ++j) {
    x_m = x_path->data[2*j];
    x_p = x_path->data[2*j+2];
    double dx = x_path->data[2*j+1] - action.getWminimum(x_m,x_p);
    double curvature = action.getWcurvature(x_m,x_p);
    S += 0.5*curvature*dx*dx + 0.5*log(curvature);
  }
  return S;
}

/* Fill in fine points */
void RotorConditionedFineAction::fill_fine_points(Path* x_path) const {
  unsigned int M_lat = x_path->M_lat;
  // interior points
  for (unsigned int j=0; j<M_lat/2-1; ++j) {
    double x_m = x_path->data[2*j];
    double x_p = x_path->data[2*(j+1)];
    double x0 = action.getWminimum(x_m,x_p);
    double sigma = 2.*action.getWcurvature(x_m,x_p);
    ExpSin2Distribution exp_sin2_dist(sigma);
    x_path->data[2*j+1] = mod_2pi(x0 + exp_sin2_dist(engine));
  }
  // Final point which requires wrap around
  double x_m = x_path->data[M_lat-2];
  double x_p = x_path->data[0];
  double x0 = action.getWminimum(x_m,x_p);
  double sigma = 2.*action.getWcurvature(x_m,x_p);
  ExpSin2Distribution exp_sin2_dist(sigma);
  x_path->data[M_lat-1] = mod_2pi(x0 + exp_sin2_dist(engine));
}

/* Evaluate conditioned action at fine points */
double RotorConditionedFineAction::evaluate(const Path* x_path) const {
  unsigned int M_lat = action.getM_lat();
  double x_m = x_path->data[M_lat-2];
  double x_p = x_path->data[0];
  double dx = x_path->data[M_lat-1] - action.getWminimum(x_m,x_p);
  double sigma = 2.0*action.getWcurvature(x_m,x_p);
  ExpSin2Distribution exp_sin2_dist(sigma);
  double S = -log(exp_sin2_dist.evaluate(dx));
  for (unsigned int j=0; j<M_lat/2-1; ++j) {
    x_m = x_path->data[2*j];
    x_p = x_path->data[2*j+2];
    double dx = x_path->data[2*j+1] - action.getWminimum(x_m,x_p);
    double sigma = 2.0*action.getWcurvature(x_m,x_p);
    ExpSin2Distribution exp_sin2_dist(sigma);
    S += -log(exp_sin2_dist.evaluate(dx));
  }
  return S;
}
