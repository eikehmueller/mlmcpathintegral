#include "gaussianconditionedfineaction.hh"
/** @file gaussianconditionedfineaction.cc
 * @brief Implementation of gaussianconditionedfineaction.hh
 */

/* Fill in fine points */
void GaussianConditionedFineAction::fill_fine_points(
    std::shared_ptr<SampleState> x_path) const {
  unsigned int M_lat = action->get_lattice()->getM_lat();
  // interior points
  for (unsigned int j = 0; j < M_lat / 2 - 1; ++j) {
    double x_m = x_path->data[2 * j];
    double x_p = x_path->data[2 * (j + 1)];
    double x0 = action->getWminimum(x_m, x_p);
    double sigma = 1. / sqrt(action->getWcurvature(x_m, x_p));
    x_path->data[2 * j + 1] = x0 + normal_dist(engine) * sigma;
  }
  // Final point which requires wrap around
  double x_m = x_path->data[M_lat - 2];
  double x_p = x_path->data[0];
  double x0 = action->getWminimum(x_m, x_p);
  double sigma = 1. / sqrt(action->getWcurvature(x_m, x_p));
  x_path->data[M_lat - 1] = x0 + normal_dist(engine) * sigma;
}

/* Evaluate conditioned action at fine points */
double GaussianConditionedFineAction::evaluate(
    const std::shared_ptr<SampleState> x_path) const {
  unsigned int M_lat = action->get_lattice()->getM_lat();
  double x_m = x_path->data[M_lat - 2];
  double x_p = x_path->data[0];
  double dx = x_path->data[M_lat - 1] - action->getWminimum(x_m, x_p);
  double curvature = action->getWcurvature(x_m, x_p);
  double S = 0.5 * curvature * dx * dx - 0.5 * log(curvature);
  for (unsigned int j = 0; j < M_lat / 2 - 1; ++j) {
    x_m = x_path->data[2 * j];
    x_p = x_path->data[2 * j + 2];
    double dx = x_path->data[2 * j + 1] - action->getWminimum(x_m, x_p);
    double curvature = action->getWcurvature(x_m, x_p);
    S += 0.5 * curvature * dx * dx - 0.5 * log(curvature);
  }
  return S;
}
