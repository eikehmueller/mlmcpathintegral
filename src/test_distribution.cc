#include <iostream>
#include <vector>
#include <memory>
#include <random>
#include <fstream>
#include <sstream>
#include <cmath>
#include "common/timer.hh"
#include "distribution/expsin2distribution.hh"
#include "distribution/besselproductdistribution.hh"

/** @file test_distribution.hh
 *
 * @brief Main program for testing probability distributions
 *
 * This file assesses different distributions by measuring the time per
 * sample/evaluation. Further, a number of samples and the value at
 * different points is written to a file which can later be analysed by a
 * Python script.
 */

/** @class DistribtionWrapper
 *
 * @brief Base class for wrapper around distribution
 *
 * This class provides a wrapper around a distribution, which allows
 * access to the following functionality:
 *
 * - Draw a single/multiple sample(s)
 * - Evaluate the distribution a single/multiple point(s)
 * - Write a header with the function parameters to an output stream
 */
class DistributionWrapper {
public:
    /** @brief Create new instance */
    DistributionWrapper() : uniform_dist(-M_PI,+M_PI) {}
    
    /** @brief Draw single sample
     *
     * @param[inout] engine Random number engibe to use
     */
    virtual double draw(std::mt19937_64& engine) = 0;
    
    /** @brief Draw multiple samples
     *
     * Use this method for time measurements to avoid overheads from
     * repeated calls.
     *
     * @param[inout] engine Random number engine to use
     * @param[in] n_samples Number of samples to draw
     */
    virtual void draw(std::mt19937_64& engine,
                             const unsigned int n_samples) = 0;

    /** @brief Evaluate at a single point
     *
     * @param[in] x Point at which the distribution is evaluated
     */
    virtual double evaluate(const double x) = 0;

    /** @brief Evaluate at multiple points
     *
     * Evaluates the function a several randomly chosen points.
     * Use this function for time measurements to avoid overheads from
     * repeated calls.
     *
     * @param[inout] engine Random number engine to use
     * @param[in] n_samples Number of points to evaluate
     */
    virtual void evaluate(std::mt19937_64& engine,
                          const unsigned int n_samples) = 0;

    /** @brief Write parameters to stream
     *
     * Write the name of the distribution and its parameters to an
     * output stream. This will be used when saving the distribution to a
     * file.
     *
     * @param[inout] out Output stream
     */
    virtual void write_header(std::ostream& out) = 0;
    
protected:
    /** @brief Uniform distribution to allow drawing from \f$[-\pi,\pi]\f$*/
    std::uniform_real_distribution<double> uniform_dist;
};

/** @class ExpSin2DistributionWrapper
 *
 * @brief Wrapper class for ExpSin2Distribution
 */
class ExpSin2DistributionWrapper : public DistributionWrapper {
public:
    /** @brief Create new instance
     *
     * @param[in] dist_ Distribution to wrap
     * @param[in] sigma_ Parameter sigma
     */
    ExpSin2DistributionWrapper(const ExpSin2Distribution& dist_,
                               const double sigma_) :
    dist(dist_), sigma(sigma_) {}
    
    /** @brief Draw single sample
     *
     * @param[inout] engine Random number engibe to use
     */
    virtual double draw(std::mt19937_64& engine) {
        return dist.draw(engine,sigma);
    }
    
    /** @brief Draw multiple samples
     *
     * Use this method for time measurements to avoid overheads from
     * repeated calls.
     *
     * @param[inout] engine Random number engine to use
     * @param[in] n_samples Number of samples to draw
     */
    virtual void draw(std::mt19937_64& engine,
                      const unsigned int n_samples) {
        for (unsigned int n=0;n<n_samples;++n) {
            double x = dist.draw(engine,sigma);
            (void) x;
        }
    }
    
    /** @brief Evaluate at a single point
     *
     * @param[in] x Point at which the distribution is evaluated
     */
    virtual double evaluate(const double x) {
        return dist.evaluate(x,sigma);
    }

    /** @brief Evaluate at multiple points
     *
     * Evaluates the function a several randomly chosen points.
     * Use this function for time measurements to avoid overheads from
     * repeated calls.
     *
     * @param[inout] engine Random number engine to use
     * @param[in] n_samples Number of points to evaluate
     */
    virtual void evaluate(std::mt19937_64& engine,
                          const unsigned int n_samples) {
        for (unsigned int j=0; j<n_samples; ++j) {
            double x = uniform_dist(engine);
            double y = dist.evaluate(x,sigma);
            (void) y;
        }
    }
   
    /** @brief Write parameters to stream
     *
     * Write the name of the distribution and its parameters to an
     * output stream. This will be used when saving the distribution to a
     * file.
     *
     * @param[inout] out Output stream
     */
    virtual void write_header(std::ostream& out) {
        out << "ExpSin2Distribution" << std::endl;
        out << "  sigma = " << sigma << std::endl;
    }
    
private:
    /** @brief Reference to distribution to wrap */
    const ExpSin2Distribution& dist;
    /** @brief Parameter \f$\sigma\f$*/
    const double sigma;
};

class BesselProductDistributionWrapper : public DistributionWrapper {
public:
    BesselProductDistributionWrapper(const BesselProductDistribution& dist_,
                                     const double x_p_,
                                     const double x_m_) :
    dist(dist_), x_p(x_p_), x_m(x_m_) {}
    
    /** @brief Draw single sample
     *
     * @param[inout] engine Random number engibe to use
     */
    virtual double draw(std::mt19937_64& engine) {
        return dist.draw(engine,x_p,x_m);
    }
    
    /** @brief Draw multiple samples
     *
     * Use this method for time measurements to avoid overheads from
     * repeated calls.
     *
     * @param[inout] engine Random number engine to use
     * @param[in] n_samples Number of samples to draw
     */
    virtual void draw(std::mt19937_64& engine,
                             const unsigned int n_samples) {
        for (unsigned int n=0;n<n_samples;++n) {
            double x = dist.draw(engine,x_p,x_m);
            (void) x;
        }
    }

    /** @brief Evaluate at a single point
     *
     * @param[in] x Point at which the distribution is evaluated
     */
    virtual double evaluate(const double x) {
        return dist.evaluate(x,x_p,x_m);
    }
    
    /** @brief Evaluate at multiple points
     *
     * Evaluates the function a several randomly chosen points.
     * Use this function for time measurements to avoid overheads from
     * repeated calls.
     *
     * @param[inout] engine Random number engine to use
     * @param[in] n_samples Number of points to evaluate
     */
    virtual void evaluate(std::mt19937_64& engine,
                                 const unsigned int n_samples) {
        for (unsigned int j=0; j<n_samples; ++j) {
            double x = uniform_dist(engine);
            double y = dist.evaluate(x,x_p,x_m);
            (void) y;
        }
    }

    /** @brief Write parameters to stream
     *
     * Write the name of the distribution and its parameters to an
     * output stream. This will be used when saving the distribution to a
     * file.
     *
     * @param[inout] out Output stream
     */
    virtual void write_header(std::ostream& out) {
        out << "BesselProductDistribution" << std::endl;
        out << "  beta = " << dist.get_beta() << std::endl;
        out << "  x_p  = " << x_p << std::endl;
        out << "  x_m  = " << x_m << std::endl;
    }
    
private:
    /** @brief Distribution to wrap */
    const BesselProductDistribution& dist;
    /** @brief Parameter \f$x_+\f$*/
    const double x_p;
    /** @brief Parameter \f$x_-\f$*/
    const double x_m;
};

/** @brief Measure time for creating a single sample
 *
 * Samples the distribution repeatedly to measure the time per sample
 *
 * @param[in] wrapper Distribution wrapper
 */
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

/** @brief Measure time for evaluating distribution at a single point
 *
 * Evaluates the distribution repeatedly to measure the time per evaluation
 *
 * @param[in] wrapper Distribution wrapper
 */
double time_evaluation(DistributionWrapper& wrapper) {
    unsigned int n_samples = 1000000;
    std::mt19937_64 engine;
    engine.seed(241857);
    std::uniform_real_distribution<double> uniform_dist(-M_PI,+M_PI);
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

/** @brief Write samples and plot of distribution to file
 *
 * This method writes a given number of samples from the distribution and
 * the value of the distribution to a file which can later be visualised with
 * a Python script. The generated file has the following format (assuming
 * n_samples = 1000, n_intervals = 128 in this example):
 *
 * -------------------------------------------------
 * BesselProductDistribution
 *   beta = 4
 *   x_p  = 2.82743
 *   x_m  = 0
 *
 * n_samples = 1000
 * n_points = 129
 *
 * ==== samples ====
 * 2.05174
 * 1.82559
 * 0.806089
 * [...]
 * 1.50355
 *
 * ==== points ====
 * -3.14159 0.0371872
 * -3.09251 0.0364034
 * [...]
 * 3.09251 0.0385538
 * 3.14159 0.0371872
 * -------------------------------------------------
 *
 * The header will vary depending on the distribution that is analysed.
 * The number of points is the number of intervals plus 1.
 *
 * @param[in] wrapper Distribution wrapper
 * @param[in] n_samples Number of samples to generate
 * @param[in] n_intervals Number of intervals used for plotting
 * @param[in] filename Name of file to write data to
 */
void save_distribution(DistributionWrapper& wrapper,
                         const unsigned int n_samples,
                         const unsigned int n_intervals,
                         const std::string filename) {
    // Open file and write header
    std::ofstream outfile(filename);
    wrapper.write_header(outfile);
    // Generate samples
    std::mt19937_64 engine;
    engine.seed(241857);
    std::vector<double> data(n_samples);
    for (unsigned int n=0; n<n_samples; ++n) {
        double x = wrapper.draw(engine);
        data[n] = x;
    }
    outfile << std::endl;
    outfile << "  n_samples = " << n_samples << std::endl;
    outfile << "  n_points = " << (n_intervals+1) << std::endl << std::endl;
    outfile << "==== samples ====" << std::endl;
    for (int n=0; n<n_samples; ++n) {
        outfile << data[n] << std::endl;
    }
    outfile << std::endl;
    outfile << "==== points ====" << std::endl;
    // Generate points
    for (int n=0;n<n_intervals+1;++n) {
        double x = -M_PI+(2.*M_PI/n_intervals)*n;
        double y = wrapper.evaluate(x);
        outfile << x << " " << y << std::endl;
    }
    outfile.close();
}

/** @brief Assess a given distribution
 *
 * Measures the time per sample/evaluation and saves results of sampling and
 * evaluation to a file which can later be analysed with a Python script.
 *
 * @param[in] wrapper Distribution wrapper
 * @param[in] n_samples Number of samples to generate
 * @param[in] n_intervals Number of intervals used for plotting
 * @param[in] filename Name of file to write data to
 */
void assess_distribution(DistributionWrapper& wrapper,
                         const unsigned int n_samples,
                         const unsigned int n_intervals,
                         const std::string filename) {
    double time_per_sample = time_sample(wrapper);
    double time_per_evaluation = time_evaluation(wrapper);
    std::cout << "time per sample = " << 1.E9*time_per_sample << " ns" << std::endl;
    std::cout << "time per evaluation = " << 1.E9*time_per_evaluation << " ns" << std::endl;
    save_distribution(wrapper,n_samples,n_intervals,filename);
}

/* *************************** M A I N ***************************** */
int main(int argc, char* argv[]) {
    // Number of samples
    unsigned int n_samples = 1000000;
    // Number of intervals for plotting the distribution
    unsigned int n_intervals = 128;
    
    /* === ExpSin2Distribution === */
    std::cout << "Testing ExpSin2Distribution ..." << std::endl;
    double sigma = 16.0;
    std::cout << "sigma = " << sigma << std::endl;
    ExpSin2Distribution expsin2_dist;
    /* Lambda expression for drawing from distribution */
    ExpSin2DistributionWrapper expsin2_wrapper(expsin2_dist,sigma);
    assess_distribution(expsin2_wrapper,
                        n_samples,
                        n_intervals,
                        "distribution_expsin2.txt");
    /* === ExpSin2Distribution === */
    std::cout << std::endl;
    
    std::cout << "Testing BesselProductDistribution ..." << std::endl;
    double beta = 4.0;
    double x_p = 0.9*M_PI;
    double x_m = 0;
    std::cout << "beta = " << beta << std::endl;
    std::cout << "x_p  = " << x_p << std::endl;
    std::cout << "x_m  = " << x_m << std::endl;
    BesselProductDistribution besselproduct_dist(beta);
    /* Lambda expression for drawing from distribution */
    BesselProductDistributionWrapper besselproduct_wrapper(besselproduct_dist,x_p,x_m);
    assess_distribution(besselproduct_wrapper,
                        n_samples,
                        n_intervals,
                        "distribution_besselproduct.txt");
}
