#include "harmonicoscillatoraction.hh"
#include <iostream>
/** @file harmonicoscillatoraction.cc
 * @brief Implementation of harmonicoscillatoraction.hh
 */

/** Evaluate action */
const double HarmonicOscillatorAction::evaluate(const Path* x_path) const {
  double ainv2 = 1./(a_lat*a_lat);
  double x_diff = x_path->data[0]-x_path->data[Action::M_lat-1];
  double S = ainv2*x_diff*x_diff + mu2*x_path->data[0]*x_path->data[0];
  for (unsigned int j=1;j<Action::M_lat;++j) {
    double x_diff = x_path->data[j]-x_path->data[j-1];
    S += ainv2*x_diff*x_diff + mu2*x_path->data[j]*x_path->data[j];
  }
  return 0.5*a_lat*m0*S;
}

/** Build Cholesky factor of covariance matrix */
void HarmonicOscillatorAction::build_covariance() {
  Matrix Sigma;
  Sigma = Matrix(Action::M_lat,Action::M_lat);
  for (unsigned int i=0; i<Action::M_lat; ++i) {
    for (unsigned int j=0; j<Action::M_lat; ++j) {
      Sigma(i,j) = 0.0;
    }
  }
  double d_tmp = a_lat*m0*mu2 + 2.0*m0/a_lat;
  double c_tmp = -m0/a_lat;
  for (unsigned int i=0; i<Action::M_lat; ++i) {
    Sigma(i,i) = d_tmp;
    Sigma(i,(i+1)%Action::M_lat) = c_tmp;
    Sigma(i,(i-1)%Action::M_lat) = c_tmp;
  }
  Eigen::LLT<Matrix> llt;
  llt.compute(Sigma.inverse());
  L_cov = llt.matrixL();
}

/** Draw sample from distribution */
const void HarmonicOscillatorAction::draw(std::vector<Path*> x_path) {
  Eigen::Map<Vector> x(x_path[0]->data,Action::M_lat);
  for (unsigned int i=0; i<Action::M_lat;++i) {
    y_tmp->data[i] = normal_dist(engine);
  }  
  Eigen::Map<Vector> y(y_tmp->data,Action::M_lat);
  x = L_cov * y;
}
