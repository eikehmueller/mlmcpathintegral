#include "harmonicoscillatoraction.hh"
#include <iostream>
/** @file harmonicoscillatoraction.cc
 * @brief Implementation of harmonicoscillatoraction.hh
 */

/** Evaluate action */
const double HarmonicOscillatorAction::evaluate(
    const std::shared_ptr<SampleState> x_path) const {
  double ainv2 = 1. / (a_lat * a_lat);
  double x_diff = x_path->data[0] - x_path->data[M_lat - 1];
  double S = ainv2 * x_diff * x_diff + mu2 * x_path->data[0] * x_path->data[0];
  for (unsigned int j = 1; j < M_lat; ++j) {
    double x_diff = x_path->data[j] - x_path->data[j - 1];
    S += ainv2 * x_diff * x_diff + mu2 * x_path->data[j] * x_path->data[j];
  }
  return 0.5 * a_lat * m0 * S;
}

/** Calculate force */
void HarmonicOscillatorAction::force(
    const std::shared_ptr<SampleState> x_path,
    std::shared_ptr<SampleState> p_path) const {
  double tmp_1 = m0 / a_lat;
  double tmp_2 = 2. + a_lat * a_lat * mu2;
  p_path->data[0] = tmp_1 * (tmp_2 * x_path->data[0] - x_path->data[M_lat - 1] -
                             x_path->data[1]);
  // Interior points
  for (unsigned int j = 1; j < M_lat - 1; ++j) {
    p_path->data[j] = tmp_1 * (tmp_2 * x_path->data[j] - x_path->data[j - 1] -
                               x_path->data[j + 1]);
  }
  p_path->data[M_lat - 1] = tmp_1 * (tmp_2 * x_path->data[M_lat - 1] -
                                     x_path->data[M_lat - 2] - x_path->data[0]);
}

/** Build Cholesky factor of covariance matrix */
void HarmonicOscillatorAction::build_covariance() {
  Matrix Sigma;
  Sigma = Matrix(M_lat, M_lat);
  for (unsigned int i = 0; i < M_lat; ++i) {
    for (unsigned int j = 0; j < M_lat; ++j) {
      Sigma(i, j) = 0.0;
    }
  }
  double d_tmp = a_lat * m0 * mu2 + 2.0 * m0 / a_lat;
  double c_tmp = -m0 / a_lat;
  for (unsigned int i = 0; i < M_lat; ++i) {
    Sigma(i, i) = d_tmp;
    Sigma(i, (i + 1) % M_lat) = c_tmp;
    Sigma(i, (i - 1) % M_lat) = c_tmp;
  }
  Eigen::LLT<Matrix> llt;
  llt.compute(Sigma.inverse());
  L_cov = llt.matrixL();
}

/** Draw sample from distribution */
void HarmonicOscillatorAction::draw(std::shared_ptr<SampleState> x_path) {
  std::generate(y_tmp->data.data(), y_tmp->data.data() + y_tmp->data.size(),
                [this]() { return normal_dist(engine); });
  x_path->data = L_cov * y_tmp->data;
  n_total_samples++;
  n_accepted_samples++;
  accept = true;
}

/** Exact expression for expectation value of \f$X^2\f$ */
const double HarmonicOscillatorAction::Xsquared_analytical() {
  double R = 1. + 0.5 * a_lat * a_lat * mu2 -
             a_lat * sqrt(mu2) * sqrt(1. + 0.25 * a_lat * a_lat * mu2);
  return 1. / (2. * m0 * sqrt(mu2) * sqrt(1 + 0.25 * a_lat * a_lat * mu2)) *
         (1. + pow(R, M_lat)) / (1. - pow(R, M_lat));
}

/** Continuum limit of expectation value of \f$X^2\f$ */
const double HarmonicOscillatorAction::Xsquared_analytical_continuum() {
  return 1. / (2. * m0 * sqrt(mu2)) * (1. + exp(-sqrt(mu2) * T_final)) /
         (1. - exp(-sqrt(mu2) * T_final));
}
