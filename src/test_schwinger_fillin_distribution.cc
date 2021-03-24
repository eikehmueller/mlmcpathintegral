#include <iostream>
#include <vector>
#include <memory>
#include <random>
#include <fstream>
#include <sstream>
#include <cmath>
#include "common/timer.hh"
#include "common/commandlineparser.hh"
#include "distribution/expcosdistribution.hh"
#include "distribution/besselproductdistribution.hh"

/** @file test_schwinger_fillin_distribution.hh
 *
 * @brief Main program for testing the fill-in distribution for the quenched  Schwinger model
 *
 * Draw samples from the distribution
 * 
 * \f[
 *  \pi(\theta_1,\theta_2,\theta_3,\theta_4)
 *  = \mathcal{N}^{-1} \exp\left[ \beta \left( \cos(\theta_1-\theta_2-\phi_{12})
 *                                           + \cos(\theta_2-\theta_3-\phi_{23})  
 *                                           + \cos(\theta_3-\theta_4-\phi_{34})  
 *                                           + \cos(\theta_4-\theta_1-\phi_{41})  
 *                                    \right) 
 *                        \right]
 * \f]
 *
 * by using the product approach and write them to a file for plotting.
 * 
 */
 
 /** @brief Generate samples and save to disk
  * 
  * After a header of the form "beta = ...", the file will contain comma
  * separated values in the following format:
  * 
  * theta_{1}^{(1)}, \theta{2}^{(1)}, \theta_{3}^{(1)}, \theta_{4}^{(1)}
  * theta_{1}^{(2)}, \theta{2}^{(2)}, \theta_{3}^{(2)}, \theta_{4}^{(2)}
  * theta_{1}^{(3)}, \theta{2}^{(3)}, \theta_{3}^{(3)}, \theta_{4}^{(3)}
  * ...
  * 
  * @param[in] beta Coupling constant \f$beta\f$
  * @param[in] n_samples Number of samples to generate
  * @param[in] phi_12 Perimeter angle \f$\phi_{12}\f$
  * @param[in] phi_23 Perimeter angle \f$\phi_{23}\f$
  * @param[in] phi_34 Perimeter angle \f$\phi_{34}\f$
  * @param[in] phi_41 Perimeter angle \f$\phi_{41}\f$
  */
 void generate_samples(const double beta,
                       const unsigned int n_samples,
                       const double phi_12,
                       const double phi_23,
                       const double phi_34,
                       const double phi_41) {
    
    std::mt19937_64 engine;
    engine.seed(215517);
    std::uniform_real_distribution<double> uniform_dist(-M_PI,M_PI);
    BesselProductDistribution besselproduct_dist(beta);
    ExpCosDistribution expcos_dist(beta);
    std::vector<double> theta(4*n_samples);
    for (unsigned int k=0;k<n_samples;++k) {
        double theta_m = -(phi_23+phi_34); 
        double theta_p = phi_12+phi_41;
        double theta_tilde = besselproduct_dist.draw(engine,theta_p,theta_m);
        double xi = uniform_dist(engine);
        double theta_2 = mod_2pi(-0.5*theta_tilde+xi);
        double theta_4 = mod_2pi(+0.5*theta_tilde+xi);
        theta_m = -phi_41+theta_4;
        theta_p = phi_12+theta_2;
        double theta_1 = expcos_dist.draw(engine,theta_p,theta_m);
        theta_m = phi_34+theta_4;
        theta_p = -phi_23+theta_2;
        double theta_3 = expcos_dist.draw(engine,theta_p,theta_m);
        theta[4*k+0] = theta_1;
        theta[4*k+1] = theta_2;
        theta[4*k+2] = theta_3;
        theta[4*k+3] = theta_4;
    }
    FILE* datafile;
    datafile = fopen("fillin_distribution.txt","w");
    fprintf(datafile,"beta = %f\n",beta);
    fprintf(datafile,"phi12 = %f\n",phi_12);
    fprintf(datafile,"phi23 = %f\n",phi_23);
    fprintf(datafile,"phi34 = %f\n",phi_34);
    fprintf(datafile,"phi41 = %f\n",phi_41);
    fprintf(datafile,"theta1,theta2,theta3,theta4\n");
    for (unsigned int k=0; k<n_samples;++k) {
        fprintf(datafile,"%+12.8e,",theta[4*k+0]);
        fprintf(datafile,"%+12.8e,",theta[4*k+1]);
        fprintf(datafile,"%+12.8e,",theta[4*k+2]);
        fprintf(datafile,"%+12.8e \n",theta[4*k+3]);
    }
    fclose(datafile);
 }

/* *************************** M A I N ***************************** */
int main(int argc, char* argv[]) {
    // Number of samples
    unsigned long n_samples = 1000000;
    // Coupling parameter \beta
    double beta = 4.0;
    // Perimeter angles
    double phi_12 = 0.0;
    double phi_23 = 0.0;
    double phi_34 = 0.0;
    double phi_41 = 0.0;
    
    // Parse command line arguments
    CommandLineParser commandlineparser(argc, argv);
    commandlineparser.getopt_double("beta",beta);
    commandlineparser.getopt_ulong("samples", n_samples);
    commandlineparser.getopt_double("phi12",phi_12);
    commandlineparser.getopt_double("phi23",phi_23);
    commandlineparser.getopt_double("phi34",phi_34);
    commandlineparser.getopt_double("phi41",phi_41);
    std::cout << "beta = " << beta << std::endl;    
    std::cout << "Number of samples = " << n_samples << std::endl;
    std::cout << "phi12 = " << phi_12 << std::endl;    
    std::cout << "phi23 = " << phi_23 << std::endl;    
    std::cout << "phi34 = " << phi_34 << std::endl;    
    std::cout << "phi41 = " << phi_41 << std::endl;    
    
    generate_samples(beta,n_samples,phi_12,phi_23,phi_34,phi_41);
    
}
