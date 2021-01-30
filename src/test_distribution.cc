#include <iostream>
#include <vector>
#include <memory>
#include <random>
#include <fstream>
#include <sstream>
#include <cmath>
#include <functional>
#include "common/timer.hh"
#include "distribution/expsin2distribution.hh"
#include "distribution/besselproductdistribution.hh"

class DistributionWrapper {
public:
    DistributionWrapper() : uniform_dist(0.,1.) {}
    
    /* Draw single sample */
    virtual double draw(std::mt19937_64 engine) = 0;
    
    /* Draw n_samples samples */
    virtual void draw(std::mt19937_64 engine,
                             const unsigned int n_samples) = 0;

    /* Evaluate n_samples samples */
    virtual void evaluate(std::mt19937_64 engine,
                                 const unsigned int n_samples) = 0;

protected:
    std::uniform_real_distribution<double> uniform_dist;
};

class ExpSin2DistributionWrapper : public DistributionWrapper {
public:
    ExpSin2DistributionWrapper(const ExpSin2Distribution& dist_,
                               const double sigma_) :
    dist(dist_), sigma(sigma_) {}
    
    /* Draw single sample */
    virtual double draw(std::mt19937_64 engine) {
        return dist.draw(engine,sigma);
    }
    
    /* Draw n_samples samples */
    virtual void draw(std::mt19937_64 engine,
                             const unsigned int n_samples) {
        for (unsigned int n=0;n<n_samples;++n) {
            double x = dist.draw(engine,sigma);
            (void) x;
        }
    }

    /* Evaluate n_samples samples */
    virtual void evaluate(std::mt19937_64 engine,
                                 const unsigned int n_samples) {
        for (unsigned int j=0; j<n_samples; ++j) {
            double x = uniform_dist(engine);
            double y = dist.evaluate(x,sigma);
            (void) y;
        }
    }
    
private:
    const ExpSin2Distribution& dist;
    const double sigma;
};

class BesselProductDistributionWrapper : public DistributionWrapper {
public:
    BesselProductDistributionWrapper(const BesselProductDistribution& dist_,
                                     const double x_p_,
                                     const double x_m_) :
    dist(dist_), x_p(x_p_), x_m(x_m_) {}
    
    /* Draw single sample */
    virtual double draw(std::mt19937_64 engine) {
        return dist.draw(engine,x_p,x_m);
    }
    
    /* Draw n_samples samples */
    virtual void draw(std::mt19937_64 engine,
                             const unsigned int n_samples) {
        for (unsigned int n=0;n<n_samples;++n) {
            double x = dist.draw(engine,x_p,x_m);
            (void) x;
        }
    }

    /* Evaluate n_samples samples */
    virtual void evaluate(std::mt19937_64 engine,
                                 const unsigned int n_samples) {
        for (unsigned int j=0; j<n_samples; ++j) {
            double x = uniform_dist(engine);
            double y = dist.evaluate(x,x_p,x_m);
            (void) y;
        }
    }
    
private:
    const BesselProductDistribution& dist;
    const double x_p;
    const double x_m;
};

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

/* Measure time for creating a single sample */
double time_sample(DistributionWrapper& wrapper) {
    unsigned int n_samples = 1000000;
    std::mt19937_64 engine;
    engine.seed(241857);
    Timer timer;
    timer.start();
    wrapper.draw(engine,n_samples);
    timer.stop();
    return timer.elapsed()/n_samples;
}

/* Measure time for evaluating distribution at a single point */
double time_evaluation(DistributionWrapper& wrapper) {
    unsigned int n_samples = 1000000;
    std::mt19937_64 engine;
    engine.seed(241857);
    std::uniform_real_distribution<double> uniform_dist(0.0,1.0);
    Timer timer;
    timer.start();
    wrapper.evaluate(engine, n_samples);
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
void binned_distribution(DistributionWrapper& wrapper,
                         const unsigned int n_samples,
                         const unsigned int n_bins,
                         const std::string filename) {
    std::mt19937_64 engine;
    engine.seed(241857);
    double rho = n_bins/(2.*M_PI*n_samples);
    std::vector<double> bin(n_bins,0);
    for (unsigned int n=0; n<n_samples; ++n) {
        double x = wrapper.draw(engine);
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

/* Evaluate a given distribution */
void evaluate_distribution(DistributionWrapper& wrapper,
                       const std::string filename) {
    double time_per_sample = time_sample(wrapper);
    double time_per_evaluation = time_evaluation(wrapper);
    std::cout << "time per sample = " << 1.E9*time_per_sample << " ns" << std::endl;
    std::cout << "time per evaluation = " << 1.E9*time_per_evaluation << " ns" << std::endl;
    unsigned int n_samples = 10000000;
    unsigned int n_bins = 32;
    binned_distribution(wrapper,n_samples,n_bins,filename);
}

/* *************************** M A I N ***************************** */
int main(int argc, char* argv[]) {
    /* === ExpSin2Distribution === */
    std::cout << "Testing ExpSin2Distribution ..." << std::endl;
    double sigma = 16.0;
    std::cout << "sigma = " << sigma << std::endl;
    ExpSin2Distribution expsin2_dist;
    /* Lambda expression for drawing from distribution */
    ExpSin2DistributionWrapper expsin2_wrapper(expsin2_dist,sigma);
    evaluate_distribution(expsin2_wrapper,"distribution_expsin2.txt");
    /* === ExpSin2Distribution === */
    std::cout << std::endl;
    std::cout << "Testing BesselProductDistribution ..." << std::endl;
    double beta = 2.0;
    double x_p = 0.2;
    double x_m = 0.7;
    std::cout << "beta = " << beta << std::endl;
    std::cout << "x_p  = " << x_p << std::endl;
    std::cout << "x_m  = " << x_m << std::endl;
    BesselProductDistribution besselproduct_dist(beta);
    /* Lambda expression for drawing from distribution */
    BesselProductDistributionWrapper besselproduct_wrapper(besselproduct_dist,x_p,x_m);
    evaluate_distribution(besselproduct_wrapper,"distribution_besselproduct.txt");
}
