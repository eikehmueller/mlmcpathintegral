#include "quenchedschwingerconditionedfineaction.hh"
/** @file quenchedschwingerconditionedfineaction.cc
 * @brief Implementation of quenchedschwingerconditionedfineaction.hh
 */

/* Fill in fine links */
void QuenchedSchwingerConditionedFineAction::fill_fine_points(std::shared_ptr<SampleState> phi_state) const {
    // Call the appropriate coarsening function
    switch (action->get_coarsening_type()) {
        case(CoarsenBoth):
            fill_fine_points_both(phi_state);
            break;
        case(CoarsenTemporal):
            fill_fine_points_temporal(phi_state);
            break;
        case(CoarsenSpatial):
            fill_fine_points_spatial(phi_state);
            break;
        default:
            mpi_parallel::cerr << "ERROR: invalid coarsening for fill_fine_points." << std::endl;
            mpi_exit(EXIT_FAILURE);
            throw std::runtime_error("...");
            break;
    }
}

/* Fill in fine links if lattice has been coarsened in both directions (2d-fill-in) */
void QuenchedSchwingerConditionedFineAction::fill_fine_points_both(std::shared_ptr<SampleState> phi_state) const {
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
    // *** Use Gaussian approximation for fill-in ***
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
    // *** Use (approximate) true distribution for fill-in ***
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
                if (not (bessel_product_dist == NULL)) {
                    // *** True distribution for vertical links ***
                    theta_tilde = bessel_product_dist->draw(engine,theta_p,theta_m);
                } else if (not (approximate_bessel_product_dist == NULL)){
                    // *** Approximate distribution for vertical links ***
                    theta_tilde = approximate_bessel_product_dist->draw(engine,theta_p,theta_m);
                } else {
                    mpi_parallel::cerr << "ERROR: no fill-in distribution defined!" << std::endl;
                    mpi_exit(EXIT_FAILURE);
                }
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

/* Fill in fine links if lattice has been coarsened in the temporal direction only (1d-fill-in) */
void QuenchedSchwingerConditionedFineAction::fill_fine_points_temporal(std::shared_ptr<SampleState> phi_state) const {
    std::shared_ptr<Lattice2D> lattice = action->get_lattice();
    const unsigned int Mt_lat = lattice->getMt_lat();
    const unsigned int Mx_lat = lattice->getMx_lat();
    /* STEP 1: Fill in temporal links on perimeter by drawing from uniform distribution */
    double dtheta;
    for (unsigned int i=0;i<Mt_lat/2;++i) {
        for (unsigned int j=0;j<Mx_lat;++j) {
            dtheta = uniform_dist(engine);
            phi_state->data[lattice->link_cart2lin(2*i  ,j  ,0)]
                = mod_2pi(phi_state->data[lattice->link_cart2lin(2*i  ,j  ,0)]+dtheta);
            phi_state->data[lattice->link_cart2lin(2*i+1,j  ,0)]
                = mod_2pi(phi_state->data[lattice->link_cart2lin(2*i+1,j  ,0)]-dtheta);
        }
    }
    /* STEP 2: Fill interior links in spatial direction */
    for (unsigned int i=0;i<Mt_lat/2;++i) {
        for (unsigned int j=0;j<Mx_lat;++j) {
            // Links pointing in temporal direction
            double theta_p  = mod_2pi(phi_state->data[lattice->link_cart2lin(2*i  ,j  ,1)]
                                    + phi_state->data[lattice->link_cart2lin(2*i  ,j+1,0)]
                                    - phi_state->data[lattice->link_cart2lin(2*i  ,j  ,0)]);
            double theta_m  = mod_2pi(phi_state->data[lattice->link_cart2lin(2*i+1,j  ,0)]
                                    + phi_state->data[lattice->link_cart2lin(2*i+2,j  ,1)]
                                    - phi_state->data[lattice->link_cart2lin(2*i+1,j+1,0)]);
            double theta_tilde = exp_cos_dist.draw(engine,theta_p,theta_m);
            phi_state->data[lattice->link_cart2lin(2*i+1,j,1)] = theta_tilde;
        }
    }

}

/* Fill in fine links if lattice has been coarsened in the spatial direction only (1d-fill-in) */
void QuenchedSchwingerConditionedFineAction::fill_fine_points_spatial(std::shared_ptr<SampleState> phi_state) const {
    std::shared_ptr<Lattice2D> lattice = action->get_lattice();
    const unsigned int Mt_lat = lattice->getMt_lat();
    const unsigned int Mx_lat = lattice->getMx_lat();
    /* STEP 1: Fill in spatial links on perimeter by drawing from uniform distribution */
    double dtheta;
    for (unsigned int i=0;i<Mt_lat;++i) {
        for (unsigned int j=0;j<Mx_lat/2;++j) {
            // Links pointing in spatial direction
            dtheta = uniform_dist(engine);
            phi_state->data[lattice->link_cart2lin(i  ,2*j  ,1)]
                = mod_2pi(phi_state->data[lattice->link_cart2lin(i  ,2*j  ,1)]+dtheta);
            phi_state->data[lattice->link_cart2lin(i  ,2*j+1,1)]
                = mod_2pi(phi_state->data[lattice->link_cart2lin(i  ,2*j+1,1)]-dtheta);
        }
    }
    /* STEP 2: Fill interior links in temporal direction */
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

/* Evaluate conditioned action at fine links */
double QuenchedSchwingerConditionedFineAction::evaluate(const std::shared_ptr<SampleState> phi_state) const {
    // Call the appropriate evaluate function
    switch (action->get_coarsening_type()) {
        case(CoarsenBoth):
            return evaluate_both(phi_state);
            break;
        case(CoarsenTemporal):
            return evaluate_temporal(phi_state);
            break;
        case(CoarsenSpatial):
            return evaluate_spatial(phi_state);
            break;
        default:
            mpi_parallel::cerr << "ERROR: invalid coarsening for evaluate." << std::endl;
            mpi_exit(EXIT_FAILURE);
            throw std::runtime_error("...");
            return 0;
            break;
    }
}


/* Evaluate conditioned action at fine links if lattice has been coarsened in both directions */
double QuenchedSchwingerConditionedFineAction::evaluate_both(const std::shared_ptr<SampleState> phi_state) const {
    std::shared_ptr<Lattice2D> lattice = action->get_lattice();
    const unsigned int Mt_lat = lattice->getMt_lat();
    const unsigned int Mx_lat = lattice->getMx_lat();
    double S = 0.0;
    if (not (gaussian_fillin_dist == NULL)) {
        // *** Use Gaussian approximation ***
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
    } else if ( not (approximate_bessel_product_dist == NULL) ) {
        // *** True distribution with approximation ***
        // Contribution from drawing vertical links
        for (unsigned int i=0;i<Mt_lat/2;++i) {
            for (unsigned int j=0;j<Mx_lat/2;++j) {
                double phi_p = mod_2pi(+ phi_state->data[lattice->link_cart2lin(2*i+1,2*j  ,0)]
                                       + phi_state->data[lattice->link_cart2lin(2*i+2,2*j  ,1)]
                                       + phi_state->data[lattice->link_cart2lin(2*i+2,2*j+1,1)]
                                       - phi_state->data[lattice->link_cart2lin(2*i+1,2*j+2,0)]);
                double phi_m = mod_2pi(- phi_state->data[lattice->link_cart2lin(2*i  ,2*j  ,0)]
                                       + phi_state->data[lattice->link_cart2lin(2*i  ,2*j  ,1)]
                                       + phi_state->data[lattice->link_cart2lin(2*i  ,2*j+1,1)]
                                       + phi_state->data[lattice->link_cart2lin(2*i  ,2*j+2,0)]);
                double theta = mod_2pi(+ phi_state->data[lattice->link_cart2lin(2*i+1,2*j  ,1)]
                                       + phi_state->data[lattice->link_cart2lin(2*i+1,2*j+1,1)]);
                S -= log(approximate_bessel_product_dist->evaluate(theta,phi_p,phi_m));
            }
        }
        // Contribution from drawing horizontal links
        for (unsigned int i=0;i<Mt_lat;++i) {
            for (unsigned int j=0;j<Mx_lat/2;++j) {
                double phi_p = mod_2pi(- phi_state->data[lattice->link_cart2lin(i  ,2*j  ,1)]
                                       + phi_state->data[lattice->link_cart2lin(i,  2*j  ,0)]
                                       + phi_state->data[lattice->link_cart2lin(i+1,2*j  ,1)]);
                double phi_m = mod_2pi(+ phi_state->data[lattice->link_cart2lin(i  ,2*j+1,1)]
                                       + phi_state->data[lattice->link_cart2lin(i,  2*j+2,0)]
                                       - phi_state->data[lattice->link_cart2lin(i+1,2*j+1,1)]);
                double theta = mod_2pi(+ phi_state->data[lattice->link_cart2lin(i  ,2*j+1,0)]);
                S -= log(exp_cos_dist.evaluate(theta,phi_p,phi_m));
            }   
        }
    } else if ( not (bessel_product_dist == NULL) ) {
        // *** True distribution ***
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
    } else {
        mpi_parallel::cerr << "ERROR: no fill-in distribution defined!" << std::endl;
        mpi_exit(EXIT_FAILURE);
    }
    return S;
}

/* Evaluate conditioned action at fine links if lattice has been coarsened in temporal direction */
double QuenchedSchwingerConditionedFineAction::evaluate_temporal(const std::shared_ptr<SampleState> phi_state) const {
    std::shared_ptr<Lattice2D> lattice = action->get_lattice();
    const unsigned int Mt_lat = lattice->getMt_lat();
    const unsigned int Mx_lat = lattice->getMx_lat();
    double S = 0.0;
    // Contribution from drawing temporal links
    for (unsigned int i=0;i<Mt_lat/2;++i) {
        for (unsigned int j=0;j<Mx_lat;++j) {
            double phi_p = mod_2pi(- phi_state->data[lattice->link_cart2lin(2*i  , j  , 0)]
                                   + phi_state->data[lattice->link_cart2lin(2*i  , j  , 1)]
                                   + phi_state->data[lattice->link_cart2lin(2*i  , j+1, 0)]);
            double phi_m = mod_2pi(+ phi_state->data[lattice->link_cart2lin(2*i+1, j  , 0)]
                                   + phi_state->data[lattice->link_cart2lin(2*i+2, j  , 1)]
                                   - phi_state->data[lattice->link_cart2lin(2*i+1, j+1, 0)]);
            double theta = mod_2pi(+ phi_state->data[lattice->link_cart2lin(2*i+1, j  , 1)]);
            S -= log(exp_cos_dist.evaluate(theta,phi_p,phi_m));
        }
    }
    return S;
}

/* Evaluate conditioned action at fine links if lattice has been coarsened in spatial direction */
double QuenchedSchwingerConditionedFineAction::evaluate_spatial(const std::shared_ptr<SampleState> phi_state) const {
    std::shared_ptr<Lattice2D> lattice = action->get_lattice();
    const unsigned int Mt_lat = lattice->getMt_lat();
    const unsigned int Mx_lat = lattice->getMx_lat();
    double S = 0.0;
    // Contribution from drawing temporal links
    for (unsigned int i=0;i<Mt_lat;++i) {
        for (unsigned int j=0;j<Mx_lat/2;++j) {
            double phi_p = mod_2pi(- phi_state->data[lattice->link_cart2lin(i  ,2*j  ,1)]
                                   + phi_state->data[lattice->link_cart2lin(i,  2*j  ,0)]
                                   + phi_state->data[lattice->link_cart2lin(i+1,2*j  ,1)]);
            double phi_m = mod_2pi(+ phi_state->data[lattice->link_cart2lin(i  ,2*j+1,1)]
                                   + phi_state->data[lattice->link_cart2lin(i,  2*j+2,0)]
                                   - phi_state->data[lattice->link_cart2lin(i+1,2*j+1,1)]);
            double theta = mod_2pi(+ phi_state->data[lattice->link_cart2lin(i  ,2*j+1,0)]);
            S -= log(exp_cos_dist.evaluate(theta,phi_p,phi_m));
        }
    }
    return S;
}
