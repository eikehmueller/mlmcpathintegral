#include "quenchedschwingerconditionedfineaction.hh"
/** @file quenchedschwingerconditionedfineaction.cc
 * @brief Implementation of quenchedschwingerconditionedfineaction.hh
 */

/* Fill in fine links */
void QuenchedSchwingerConditionedFineAction::fill_fine_points(const std::shared_ptr<SampleState> phi_state_n,
                                                              std::shared_ptr<SampleState> phi_state) const {
    std::shared_ptr<Lattice2D> lattice = action->get_lattice();
    const unsigned int Mt_lat = lattice->getMt_lat();
    const unsigned int Mx_lat = lattice->getMx_lat();
    /* STEP 1: Fill in links on perimeter by drawing from uniform distribution */
    double dtheta;
    for (unsigned int i=0;i<Mt_lat/2;++i) {
        for (unsigned int j=0;j<Mx_lat/2;++j) {
            // Links pointing in temporal direction
            dtheta = uniform_dist(engine);
            phi_state->data[lattice->link_cart2lin(2*i  ,2*j  ,0)]
                = mod_2pi(phi_state->data[lattice->link_cart2lin(2*i  ,2*j  ,0)]+dtheta);
            phi_state->data[lattice->link_cart2lin(2*i+1,2*j  ,0)]
                = mod_2pi(phi_state->data[lattice->link_cart2lin(2*i+1,2*j  ,0)]-dtheta);
            // Links pointing in spatial direction
            dtheta = uniform_dist(engine);
            phi_state->data[lattice->link_cart2lin(2*i  ,2*j  ,1)]
                = mod_2pi(phi_state->data[lattice->link_cart2lin(2*i  ,2*j  ,1)]+dtheta);
            phi_state->data[lattice->link_cart2lin(2*i  ,2*j+1,1)]
                = mod_2pi(phi_state->data[lattice->link_cart2lin(2*i  ,2*j+1,1)]-dtheta);
        }
    }
    if (not (gaussian_fillin_dist == NULL)) {
        /* STEP 2 & 3: fill in all interior links with Gaussian approximation */
        for (unsigned int i=0;i<Mt_lat/2;++i) {
            for (unsigned int j=0;j<Mx_lat/2;++j) {
                double phi_12 = mod_2pi(+phi_state->data[lattice->link_cart2lin(2*i  ,2*j+1,1)]
                                        +phi_state->data[lattice->link_cart2lin(2*i  ,2*j+2,0)]);
                double phi_23 = mod_2pi(+phi_state->data[lattice->link_cart2lin(2*i+1,2*j+2,0)]
                                        -phi_state->data[lattice->link_cart2lin(2*i+2,2*j+1,1)]);
                double phi_34 = mod_2pi(-phi_state->data[lattice->link_cart2lin(2*i+2,2*j  ,1)]
                                        -phi_state->data[lattice->link_cart2lin(2*i+1,2*j  ,0)]);
                double phi_41 = mod_2pi(-phi_state->data[lattice->link_cart2lin(2*i  ,2*j  ,0)]
                                        +phi_state->data[lattice->link_cart2lin(2*i  ,2*j  ,1)]);
                /*
                double phi_n_12 = mod_2pi(+phi_state_n->data[lattice->link_cart2lin(2*i  ,2*j+1,1)]
                                          +phi_state_n->data[lattice->link_cart2lin(2*i  ,2*j+2,0)]);
                double phi_n_23 = mod_2pi(+phi_state_n->data[lattice->link_cart2lin(2*i+1,2*j+2,0)]
                                          -phi_state_n->data[lattice->link_cart2lin(2*i+2,2*j+1,1)]);
                double phi_n_34 = mod_2pi(-phi_state_n->data[lattice->link_cart2lin(2*i+2,2*j  ,1)]
                                          -phi_state_n->data[lattice->link_cart2lin(2*i+1,2*j  ,0)]);
                double phi_n_41 = mod_2pi(-phi_state_n->data[lattice->link_cart2lin(2*i  ,2*j  ,0)]
                                          +phi_state_n->data[lattice->link_cart2lin(2*i  ,2*j  ,1)]);
                double theta_n_1 = +phi_state_n->data[lattice->link_cart2lin(2*i  ,2*j+1,0)];
                double theta_n_2 = -phi_state_n->data[lattice->link_cart2lin(2*i+1,2*j+1,1)];
                double theta_n_3 = -phi_state_n->data[lattice->link_cart2lin(2*i+1,2*j+1,0)];
                double theta_n_4 = +phi_state_n->data[lattice->link_cart2lin(2*i+1,2*j  ,1)];                
                gaussian_fillin_dist->draw_pcn(engine,
                                               phi_n_12,phi_n_23,phi_n_34,phi_n_41,
                                               theta_n_1,theta_n_2,theta_n_3,theta_n_4,
                                               phi_12,phi_23,phi_34,phi_41,
                                               theta_1,theta_2,theta_3,theta_4);
                */
                double theta_1, theta_2, theta_3, theta_4;
                gaussian_fillin_dist->draw(engine,
                                           phi_12,phi_23,phi_34,phi_41,
                                           theta_1,theta_2,theta_3,theta_4);
                phi_state->data[lattice->link_cart2lin(2*i  ,2*j+1,0)] = +theta_1;
                phi_state->data[lattice->link_cart2lin(2*i+1,2*j+1,1)] = -theta_2;
                phi_state->data[lattice->link_cart2lin(2*i+1,2*j+1,0)] = -theta_3;
                phi_state->data[lattice->link_cart2lin(2*i+1,2*j  ,1)] = +theta_4;
            }
        }
    } else {
        /* STEP 2: Fill interior links in spatial direction */
        for (unsigned int i=0;i<Mt_lat/2;++i) {
            for (unsigned int j=0;j<Mx_lat/2;++j) {
                // Links pointing in temporal direction
                double theta_p = mod_2pi(phi_state->data[lattice->link_cart2lin(2*i+1,2*j  ,0)]
                                       + phi_state->data[lattice->link_cart2lin(2*i+2,2*j  ,1)]
                                       + phi_state->data[lattice->link_cart2lin(2*i+2,2*j+1,1)]
                                       - phi_state->data[lattice->link_cart2lin(2*i+1,2*j+2,0)]);
                double theta_m = mod_2pi(phi_state->data[lattice->link_cart2lin(2*i  ,2*j  ,1)]
                                       + phi_state->data[lattice->link_cart2lin(2*i  ,2*j+1,1)]
                                       + phi_state->data[lattice->link_cart2lin(2*i  ,2*j+2,0)]
                                       - phi_state->data[lattice->link_cart2lin(2*i  ,2*j  ,0)]);
                double theta_tilde;
                theta_tilde = bessel_product_dist->draw(engine,theta_p,theta_m);
                dtheta = uniform_dist(engine);
                phi_state->data[lattice->link_cart2lin(2*i+1,2*j  ,1)] = mod_2pi(0.5*theta_tilde+dtheta);
                phi_state->data[lattice->link_cart2lin(2*i+1,2*j+1,1)] = mod_2pi(0.5*theta_tilde-dtheta);
            }
        }
        /* STEP 3: Fill interior links in temporal direction */
        for (unsigned int i=0;i<Mt_lat;++i) {
            for (unsigned int j=0;j<Mx_lat/2;++j) {
                // Links pointing in temporal direction
                double theta_p  = mod_2pi(phi_state->data[lattice->link_cart2lin(i  ,2*j  ,0)]
                                        + phi_state->data[lattice->link_cart2lin(i+1,2*j  ,1)]
                                        - phi_state->data[lattice->link_cart2lin(i  ,2*j  ,1)]);
                double theta_m  = mod_2pi(phi_state->data[lattice->link_cart2lin(i  ,2*j+1,1)]
                                        + phi_state->data[lattice->link_cart2lin(i  ,2*j+2,0)]
                                        - phi_state->data[lattice->link_cart2lin(i+1,2*j+1,1)]);
                double theta_tilde = exp_cos_dist.draw(engine,theta_p,theta_m);
                phi_state->data[lattice->link_cart2lin(i,2*j+1,0)] = theta_tilde;
            }
        }
    }
}

/* Evaluate conditioned action at fine links */
double QuenchedSchwingerConditionedFineAction::evaluate(const std::shared_ptr<SampleState> phi_state) const {
    std::shared_ptr<Lattice2D> lattice = action->get_lattice();
    const unsigned int Mt_lat = lattice->getMt_lat();
    const unsigned int Mx_lat = lattice->getMx_lat();
    double S = 0.0;
    if (not (gaussian_fillin_dist == NULL)) {
        // Use Gaussian approximation
        for (unsigned int i=0;i<Mt_lat/2;++i) {
            for (unsigned int j=0;j<Mx_lat/2;++j) {
                double phi_12 = mod_2pi(+phi_state->data[lattice->link_cart2lin(2*i  ,2*j+1,1)]
                                        +phi_state->data[lattice->link_cart2lin(2*i  ,2*j+2,0)]);
                double phi_23 = mod_2pi(+phi_state->data[lattice->link_cart2lin(2*i+1,2*j+2,0)]
                                        -phi_state->data[lattice->link_cart2lin(2*i+2,2*j+1,1)]);
                double phi_34 = mod_2pi(-phi_state->data[lattice->link_cart2lin(2*i+2,2*j  ,1)]
                                        -phi_state->data[lattice->link_cart2lin(2*i+1,2*j  ,0)]);
                double phi_41 = mod_2pi(-phi_state->data[lattice->link_cart2lin(2*i  ,2*j  ,0)]
                                        +phi_state->data[lattice->link_cart2lin(2*i  ,2*j  ,1)]);
                double theta_1 = mod_2pi(+phi_state->data[lattice->link_cart2lin(2*i  ,2*j+1,0)]);
                double theta_2 = mod_2pi(-phi_state->data[lattice->link_cart2lin(2*i+1,2*j+1,1)]);
                double theta_3 = mod_2pi(-phi_state->data[lattice->link_cart2lin(2*i+1,2*j+1,0)]);
                double theta_4 = mod_2pi(+phi_state->data[lattice->link_cart2lin(2*i+1,2*j  ,1)]);
                S -= log(gaussian_fillin_dist->evaluate(theta_1,theta_2,theta_3,theta_4,
                                                        phi_12,phi_23,phi_34,phi_41));
            }
        }
    } else {
        // Use full distribution
        for (unsigned int i=0;i<Mt_lat/2;++i) {
            for (unsigned int j=0;j<Mx_lat/2;++j) {
                double phi_12 = + phi_state->data[lattice->link_cart2lin(2*i,  2*j+1,1)]
                                + phi_state->data[lattice->link_cart2lin(2*i,  2*j+2,0)];
                double phi_23 = + phi_state->data[lattice->link_cart2lin(2*i+1,2*j+2,0)]
                                - phi_state->data[lattice->link_cart2lin(2*i+2,2*j+1,1)];
                double phi_34 = - phi_state->data[lattice->link_cart2lin(2*i+1,2*j  ,0)]
                                - phi_state->data[lattice->link_cart2lin(2*i+2,2*j  ,1)];
                double phi_41 = - phi_state->data[lattice->link_cart2lin(2*i  ,2*j  ,0)]
                                + phi_state->data[lattice->link_cart2lin(2*i  ,2*j  ,1)];
                double theta_1 = +phi_state->data[lattice->link_cart2lin(2*i  ,2*j+1,0)];
                double theta_2 = -phi_state->data[lattice->link_cart2lin(2*i+1,2*j+1,1)];
                double theta_3 = -phi_state->data[lattice->link_cart2lin(2*i+1,2*j+1,0)];
                double theta_4 = +phi_state->data[lattice->link_cart2lin(2*i+1,2*j  ,1)];
                double Phi = phi_12+phi_23+phi_34+phi_41;
                S -= beta * (  cos(theta_1-theta_2-phi_12)
                             + cos(theta_2-theta_3-phi_23)
                             + cos(theta_3-theta_4-phi_34)
                             + cos(theta_4-theta_1-phi_41) );
                S -= log(bessel_product_dist->Znorm_inv(Phi,true));
            }
        }
    }
    return S;
}
