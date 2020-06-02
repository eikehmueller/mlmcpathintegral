#include "conditionedfineaction.hh"
/** @file conditionedfineaction->cc
 * @brief Implementation of conditionedfineaction->hh
 */

/* Fill in fine points */
void GaussianConditionedFineAction::fill_fine_points(std::shared_ptr<Path> x_path) const {
  unsigned int M_lat = x_path->M_lat;
  // interior points
  for (unsigned int j=0; j<M_lat/2-1; ++j) {
    double x_m = x_path->data[2*j];
    double x_p = x_path->data[2*(j+1)];
    double x0 = action->getWminimum(x_m,x_p);
    double sigma = 1./sqrt(action->getWcurvature(x_m,x_p));
    x_path->data[2*j+1] = x0 + normal_dist(engine)*sigma;
  }
  // Final point which requires wrap around
  double x_m = x_path->data[M_lat-2];
  double x_p = x_path->data[0];
  double x0 = action->getWminimum(x_m,x_p);
  double sigma = 1./sqrt(action->getWcurvature(x_m,x_p));
  x_path->data[M_lat-1] = x0 + normal_dist(engine)*sigma;
}

/* Evaluate conditioned action at fine points */
double GaussianConditionedFineAction::evaluate(const std::shared_ptr<Path> x_path) const {
  unsigned int M_lat = action->getM_lat();
  double x_m = x_path->data[M_lat-2];
  double x_p = x_path->data[0];
  double dx = x_path->data[M_lat-1] - action->getWminimum(x_m,x_p);
  double curvature = action->getWcurvature(x_m,x_p);
  double S = 0.5*curvature*dx*dx - 0.5*log(curvature);
  for (unsigned int j=0; j<M_lat/2-1; ++j) {
    x_m = x_path->data[2*j];
    x_p = x_path->data[2*j+2];
    double dx = x_path->data[2*j+1] - action->getWminimum(x_m,x_p);
    double curvature = action->getWcurvature(x_m,x_p);
    S += 0.5*curvature*dx*dx - 0.5*log(curvature);
  }
  return S;
}

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

/* Find midpoint, taking into account periodicity */
double RotorNonTrigonometricConditionedFineAction::midpoint(const double x_m, const double x_p) const {
  double delta[3];
  double mp[3];
  delta[0] = fabs(x_m-x_p);
  delta[1] = fabs(x_p-x_m+2.*M_PI);
  delta[2] = fabs(x_p-x_m-2.*M_PI);
  mp[0] = 0.5*(x_m+x_p);
  mp[1] = mp[0]+M_PI;
  mp[2] = mp[0]-M_PI;
  return mod_2pi(mp[std::distance(delta,std::min_element(delta,delta+3))]);
}

/* Fill in fine points */
void RotorNonTrigonometricConditionedFineAction::fill_fine_points(std::shared_ptr<Path> x_path) const {
  unsigned int M_lat = x_path->M_lat;
  double sigma = sqrt(0.5*action->geta_lat()/action->getm0());
  // interior points
  for (unsigned int j=0; j<M_lat/2-1; ++j) {
    double x_m = x_path->data[2*j];
    double x_p = x_path->data[2*(j+1)];
    double x0 = midpoint(x_m,x_p);
    x_path->data[2*j+1] = mod_2pi(x0 + sigma*norm_dist(engine));
  }
  // Final point which requires wrap around
  double x_m = x_path->data[M_lat-2];
  double x_p = x_path->data[0];
  double x0 = midpoint(x_m,x_p);
  x_path->data[M_lat-1] = mod_2pi(x0 + sigma*norm_dist(engine));
}

/* Evaluate conditioned action at fine points */
double RotorNonTrigonometricConditionedFineAction::evaluate(const std::shared_ptr<Path> x_path) const {
  unsigned int M_lat = action->getM_lat();
  int k_max=8;
  double twosigma2_inv = action->getm0()/action->geta_lat();
  double x_m = x_path->data[M_lat-2];
  double x_p = x_path->data[0];
  double dx = x_path->data[M_lat-1] - midpoint(x_m,x_p);
  double p_local = 0.0;
  for (int k=-k_max;k<=k_max;++k) {
    double dx_tmp = dx-2.0*M_PI*k;
    p_local += exp(-twosigma2_inv*dx_tmp*dx_tmp);
  }
  double S = -log(p_local);
  for (unsigned int j=0; j<M_lat/2-1; ++j) {
    x_m = x_path->data[2*j];
    x_p = x_path->data[2*j+2];
    double dx = x_path->data[2*j+1] - midpoint(x_m,x_p);
    p_local = 0.0;
    for (int k=-k_max;k<=k_max;++k) {
      double dx_tmp = dx-2.0*M_PI*k;
      p_local += exp(-twosigma2_inv*dx_tmp*dx_tmp);
    }
    S += -log(p_local);
  }
  return S;
}
