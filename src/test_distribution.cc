#include <iostream>
#include <vector>
#include <memory>
#include <random>
#include <fstream>
#include <sstream>
#include <cmath>
#include "common/timer.hh"
#include "distribution/expsin2distribution.hh"

/* Measure time for creating a single sample */
double time_sample(const ExpSin2Distribution& dist,
                   const double sigma) {
    unsigned int n_samples = 1000000;
    std::mt19937_64 engine;
    engine.seed(241857);
    Timer timer;
    timer.start();
    for (unsigned int j=0; j<n_samples; ++j) {
        double x = dist.draw(engine,sigma);
        (void) x;
    }
    timer.stop();
    return timer.elapsed()/n_samples;
}

/* Measure time for evaluating distribution at a single point */
double time_evaluation(const ExpSin2Distribution& dist,
                       const double sigma) {
    unsigned int n_samples = 1000000;
    std::mt19937_64 engine;
    engine.seed(241857);
    std::uniform_real_distribution<double> uniform_dist(0.0,1.0);
    Timer timer;
    timer.start();
    for (unsigned int j=0; j<n_samples; ++j) {
        double x = uniform_dist(engine);
        double y = dist.evaluate(x,sigma);
        (void) y;
    }
    timer.stop();
    Timer timer_uniform_dist;
    timer_uniform_dist.start();
    for (unsigned int j=0; j<n_samples; ++j) {
        double x = uniform_dist(engine);
        (void) x;
    }
    timer_uniform_dist.stop();
    return (timer.elapsed()-timer_uniform_dist.elapsed())/n_samples;
}

/* Create binned distribution */
void binned_distribution(const ExpSin2Distribution& dist,
                         const double sigma,
                         const unsigned int n_samples,
                         const unsigned int n_bins,
                         const std::string filename) {
    std::mt19937_64 engine;
    engine.seed(241857);
    double rho = n_bins/(2.*M_PI*n_samples);
    std::vector<double> bin(n_bins,0);
    for (unsigned int n=0; n<n_samples; ++n) {
        double x = dist.draw(engine,sigma);
        int k = int(floor(n_bins*(x+M_PI)/(2.*M_PI)));
        bin[k]+=rho;
    }
    std::ofstream outfile(filename);
    for (int k=0; k<n_bins; ++k) {
        outfile << bin[k] << " ";
    }
    outfile << std::endl;
    outfile.close();
}

/* *************************** M A I N ***************************** */
int main(int argc, char* argv[]) {
    std::cout << "Testing ExpSin2Distribution ..." << std::endl;
    double sigma = 16.0;
    std::cout << "sigma = " << sigma << std::endl;
    ExpSin2Distribution dist;
    double time_per_sample = time_sample(dist,sigma);
    double time_per_evaluation = time_evaluation(dist,sigma);
    std::cout << "time per sample = " << 1.E9*time_per_sample << " ns" << std::endl;
    std::cout << "time per evaluation = " << 1.E9*time_per_evaluation << " ns" << std::endl;
    unsigned int n_samples = 10000000;
    unsigned int n_bins = 32;
    std::string filename = "distribution.txt";
    binned_distribution(dist,sigma,n_samples,n_bins,filename);
}
