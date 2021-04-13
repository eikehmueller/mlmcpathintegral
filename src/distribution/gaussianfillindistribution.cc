#include "gaussianfillindistribution.hh"
/** @file gaussianfillindistribution.cc
 * @brief Implementation of gaussianfillindistribution.hh
 */

/* Evaluate at a given point */
double GaussianFillinDistribution::evaluate(const double theta_1, const double theta_2,
                                            const double theta_3, const double theta_4,
                                            const double phi_12, const double phi_23,
                                            const double phi_34, const double phi_41) const {
    double eta_1 = mod_2pi(0.5*(theta_1+theta_2-theta_3-theta_4) + 0.5*(phi_41-phi_23));
    double eta_2 = mod_2pi(0.5*(theta_1-theta_2-theta_3+theta_4) + 0.5*(phi_34-phi_12));
    double eta_3 = mod_2pi(0.5*(theta_1-theta_2+theta_3-theta_4) + 0.25*(-phi_12+phi_23-phi_34+phi_41));
    double Phi = 0.25*(phi_12+phi_23+phi_34+phi_41);
    bool swap_eta = false;
    bool shift_eta = false;
    double Phi_star = Phi;
    if (Phi < 0.) {
        Phi_star *= -1.0;
        swap_eta = true;
    }
    if (Phi_star > 0.5*M_PI) {
        Phi_star = M_PI-Phi_star;
        swap_eta = not(swap_eta);
        shift_eta = true;
    }
    if (swap_eta) {
        std::swap(eta_1,eta_2);        
    }
    if (shift_eta) {
        eta_1 = mod_2pi(eta_1+M_PI);
        eta_2 = mod_2pi(eta_2+M_PI);
    }
    // Main peak
    double g_c = 0.0;    
    double mu_c_1[9] = { 0, 1, 1, 1, 1,-1,-1,-1,-1};
    double mu_c_2[9] = { 0, 1, 1,-1,-1, 1, 1,-1,-1};
    double mu_c_3[9] = { 0, 1,-1, 1,-1, 1,-1, 1,-1};
    double sigma2_inv_c = 2.*beta*cos(Phi_star);
    for (int j=0;j<9;++j) {
        double d_eta_1 = eta_1 - mu_c_1[j]*M_PI;
        double d_eta_2 = eta_2 - mu_c_2[j]*M_PI;
        double d_eta_3 = eta_3 - mu_c_3[j]*M_PI;
        double Q = d_eta_1*d_eta_1 + d_eta_2*d_eta_2 + 2.*d_eta_3*d_eta_3;
        g_c += exp(-0.5*sigma2_inv_c*Q);
    }
    // Secondary peak
    double g_s = 0.0;
    double mu_s_1[4] = { 0.0, 0.0,-1.0, 1.0};
    double mu_s_2[4] = {-1.0, 1.0, 0.0, 0.0};
    double mu_s_3[4] = {-0.5,-0.5, 0.5, 0.5};
    double sigma2_inv_s = 2.*beta*sin(Phi_star);
    for (int j=0;j<4;++j) {
        double d_eta_1 = eta_1 - mu_s_1[j]*M_PI;
        double d_eta_2 = eta_2 - mu_s_2[j]*M_PI;
        double d_eta_3 = eta_3 - mu_s_3[j]*M_PI;
        double Q = d_eta_1*d_eta_1 + d_eta_2*d_eta_2 + 2.*d_eta_3*d_eta_3;
        g_s += exp(-0.5*sigma2_inv_s*Q);
    }
    double p_c = get_pc(Phi_star);
    double norm_c = pow(sigma2_inv_c,1.5);
    double norm_s = pow(sigma2_inv_s,1.5);
    return p_c*norm_c*g_c + (1.-p_c)*norm_s*g_s;
}
