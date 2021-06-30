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
    double Phi_star = Phi;
    if (Phi_star < 0.) {
        Phi_star *= -1.0;
        swap_eta = true;
    }
    if (Phi_star > 0.5*M_PI) {
        Phi_star = M_PI-Phi_star;
        swap_eta = not(swap_eta);        
        eta_1 = mod_2pi(eta_1+M_PI);
        eta_2 = mod_2pi(eta_2+M_PI);
    }
    if (swap_eta) {
        std::swap(eta_1,eta_2);        
    }
    double p_c = get_pc(Phi_star);
    // If we did not add Gaussian noise, we know that eta = (0,0,0) for the
    // main peak and eta != (0,0,0) for the secondary peak. Otherwise, we need to
    // evaluate the pdf, which is the sum of two peaks.
    if (add_gaussian_noise) {
        // Evaluate normal distributions at both peaks
        double g_c = 0.0; 
        double g_s = 0.0;   
        double sigma2_inv_c = 2.*beta*cos(Phi_star);
        double sigma2_inv_s = 2.*beta*sin(Phi_star);
        // Main peaks
        for (auto p=main_peaks.begin();p!=main_peaks.end();++p) {            
            double d_eta_1 = eta_1 - (*p)[0];
            double d_eta_2 = eta_2 - (*p)[1];
            double d_eta_3 = eta_3 - (*p)[2];
            double Q = d_eta_1*d_eta_1 + d_eta_2*d_eta_2 + 2.*d_eta_3*d_eta_3;
            g_c += exp(-0.5*sigma2_inv_c*Q);        
        }
        // Secondary peaks
        for (auto p=secondary_peaks.begin();p!=secondary_peaks.end();++p) {
            double d_eta_1 = eta_1 - (*p)[0];
            double d_eta_2 = eta_2 - (*p)[1];
            double d_eta_3 = eta_3 - (*p)[2];
            double Q = d_eta_1*d_eta_1 + d_eta_2*d_eta_2 + 2.*d_eta_3*d_eta_3;
            g_s += exp(-0.5*sigma2_inv_s*Q);
        }
        // normalisation constants
        double norm_c = pow(sigma2_inv_c,1.5);
        double norm_s = pow(sigma2_inv_s,1.5);
        return p_c*norm_c*g_c + (1.-p_c)*norm_s*g_s;
    } else {
        double tolerance = 1.E-12;
        if (eta_1*eta_1+eta_2*eta_2+eta_3*eta_3 < tolerance) {
            // Main peak
            return p_c;
        } else {
            // Secondary peak
            return 1-p_c;
        }
    }
}

/* Construct peak locations */
void GaussianFillinDistribution::construct_peaks(int n_offsets) {
    // Construct unique set of peak positions in units of \pi/2
    std::set<std::array<int,3> > main_indices;
    std::set<std::array<int,3> > secondary_indices;
    // Position of peaks in unit box
    std::array<std::array<int,3>,9> p_main = { 0, 0, 0,
                                               2, 2, 2,
                                              -2, 2, 2,
                                               2,-2, 2,
                                              -2,-2, 2,
                                               2, 2,-2,
                                              -2, 2,-2,
                                               2,-2,-2,
                                              -2,-2,-2};
    std::array<std::array<int,3>,4> p_secondary = { 2, 0, 1,
                                                   -2, 0, 1,
                                                    0, 2,-1,
                                                    0,-2,-1};
    for (int k_x=-n_offsets;k_x<=+n_offsets;++k_x) {
        for (int k_y=-n_offsets;k_y<=+n_offsets;++k_y) {
            for (int k_z=-n_offsets;k_z<=+n_offsets;++k_z) {
                // Add offset
                auto add_offset = [&](std::array<int,3> a){
                    return std::array<int,3> { a[0]+4*k_x,
                                               a[1]+4*k_y,
                                               a[2]+4*k_z}; };
                for (auto x=p_main.begin();x!=p_main.end();++x) {
                    main_indices.insert(add_offset(*x));
                }
                for (auto x=p_secondary.begin();x!=p_secondary.end();++x) {
                    secondary_indices.insert(add_offset(*x));
                }
            }
        }
    }
    // Scale by multiples of \pi/2 to convert back to units of standard box
    main_peaks.clear();
    for (auto x=main_indices.begin();x!=main_indices.end();++x) {
        std::array<double,3> a;
        for (int j=0;j<3;++j) a[j]=0.5*M_PI*(*x)[j];
        main_peaks.push_back(a);			 
    }
    secondary_peaks.clear();
    for (auto x=secondary_indices.begin();x!=secondary_indices.end();++x) {
        std::array<double,3> a;
        for (int j=0;j<3;++j) a[j]=0.5*M_PI*(*x)[j];
        secondary_peaks.push_back(a);
    }
}
