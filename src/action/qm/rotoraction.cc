#include "rotoraction.hh"

/** @file rotoraction.cc
 * @brief Implementation of rotoraction.hh
 */

/** Evaluate action */
const double RotorAction::evaluate(const std::shared_ptr<SampleState> x_path) const {
    double x_diff = x_path->data[0]-x_path->data[M_lat-1];
    double S = 1.-cos(x_diff);
    for (unsigned int j=1; j<M_lat; ++j) {
        double x_diff = x_path->data[j]-x_path->data[j-1];
        S += 1.-cos(x_diff);
    }
    return m0/a_lat*S;
}

/** Local heat bath update */
void RotorAction::heatbath_update(std::shared_ptr<SampleState> x_path,
                                  const unsigned int j) {
  double x_m;
  double x_p;
  if (j==0) {
    x_m = x_path->data[M_lat-1];
  } else {
    x_m = x_path->data[j-1];
  }
  if (j==M_lat-1) {
    x_p = x_path->data[0];
  } else {
    x_p = x_path->data[j+1];
  }
  double x0 = getWminimum(x_m,x_p);
  double sigma = 2.*getWcurvature(x_m,x_p);
  x_path->data[j] = mod_2pi(x0 + exp_sin2_dist.draw(engine,sigma));
}

/** Local overrelaxation update */
void RotorAction::overrelaxation_update(std::shared_ptr<SampleState> x_path, const unsigned int j) {
    double x_m;
    double x_p;
    if (j==0) {
      x_m = x_path->data[M_lat-1];
    } else {
      x_m = x_path->data[j-1];
    }
    if (j==M_lat-1) {
      x_p = x_path->data[0];
    } else {
      x_p = x_path->data[j+1];
    }
    double x0 = getWminimum(x_m,x_p);
    x_path->data[j] = mod_2pi(2.0*x0 - x_path->data[j]);
}

/** Calculate force */
void RotorAction::force(const std::shared_ptr<SampleState> x_path,
                        std::shared_ptr<SampleState> p_path) const {
    double tmp = m0/a_lat;
    // Left boundary
    double x_m = x_path->data[M_lat-1];
    double x = x_path->data[0];
    double x_p = x_path->data[1];
    p_path->data[0] = tmp*(sin(x-x_m)+sin(x-x_p));
    // Interior points
    for (unsigned int j=1; j<M_lat-1; ++j) {
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
void RotorAction::initialise_state(std::shared_ptr<SampleState> x_path) const {
    std::uniform_real_distribution<double> uniform(-M_PI,M_PI);
    std::generate(x_path->data.data(),
                  x_path->data.data()+x_path->data.size(),
    [this,&uniform]() {
        return uniform(engine);
    });
#ifdef SAVE_PATHS
    x_path->save_to_disk("path_initial.dat");
#endif // SAVE_PATHS
}

/** Exact analytical expression for chi_t for finite a */
double RotorAction::chit_exact() const {
    double T_final = lattice->getT_final();
    return 1./m0*Phi_chit(m0/a_lat,T_final/a_lat);
}


/** Perturbative expression for chi_t for finite a */
double RotorAction::chit_perturbative() const {
    double T_final = lattice->getT_final();
    double xi = T_final/m0;
    double z = a_lat/m0;
    double S_hat2 = Sigma_hat(xi,2);
    double S_hat4 = Sigma_hat(xi,4);
    return 1./(4.*M_PI*M_PI*m0)*(1.-xi*S_hat2+(0.5-xi*S_hat2+0.25*xi*xi*(S_hat4-S_hat2*S_hat2))*z);
}

/** Analytical expression for chi_t in continuum limit */
double RotorAction::chit_continuum() const {
    double T_final = lattice->getT_final();
    double xi = T_final/m0;
    return 1./(4.*M_PI*M_PI*m0)*(1.-xi*Sigma_hat(xi,2));
}
