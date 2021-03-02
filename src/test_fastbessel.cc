#include <stdio.h>
#include <iostream>
#include <random>
#include <gsl/gsl_sf_bessel.h>
#include "common/timer.hh"
#include "common/fastbessel.hh"

template <int N>
struct ModifiedBesselCoefficientPrinter
{
    static void show() {
        ModifiedBesselCoefficientPrinter<N-1>::show();
        printf("a[%3d] = ",N);
        printf("%14.12e",ModifiedBesselCoefficient<N>::value);
        std::cout << std::endl;
    }
};

template <>
struct ModifiedBesselCoefficientPrinter<-1>
{
    static void show() {};
};

double time_evaluation(const double z_min,
                       const double z_max) {
    unsigned long n_samples = 50000000;
    std::mt19937_64 engine;
    engine.seed(241857);
    std::uniform_real_distribution<double> uniform_dist(z_min,z_max);
    Timer timer;
    timer.start();
    for (unsigned int j=0; j<n_samples; ++j) {
        double z = uniform_dist(engine);
        double x = fast_bessel_I0_scaled(z);
        (void) x;
    }
    timer.stop();
    Timer timer_uniform_dist;
    timer_uniform_dist.start();
    for (unsigned int j=0; j<n_samples; ++j) {
        double z = uniform_dist(engine);
        (void) z;
    }
    timer_uniform_dist.stop();
    return (timer.elapsed()-timer_uniform_dist.elapsed())/n_samples;
}

int main(int argc, char* argv[]) {
    std::cout << "*** expansion coefficients ***" << std::endl;
    ModifiedBesselCoefficientPrinter<16>::show();
    FILE * file_id;
    file_id = fopen ("fastbessel.csv","w");
    for (int k=0;k<2000;++k) {
        double z =k;
        double i0_scaled_exact = gsl_sf_bessel_I0_scaled(z);
        double i0_scaled_fast = fast_bessel_I0_scaled(z);
        double error = fabs(i0_scaled_fast-i0_scaled_exact)/i0_scaled_exact;
        fprintf(file_id,"%12.4f, %8.4e\n",z,error);
    }
    fclose(file_id);
    double range[5][2] = {  {   1.0, 100.0},
                            { 101.0, 200.0},
                            { 201.0, 400.0},
                            { 401.0,1100.0},
                            {1101.0,2000.0} };
    std::cout << std::endl;
    std::cout << "*** time per evaluation ***" << std::endl;
    for (size_t k=0;k<5;++k) {
        double z_min = range[k][0];
        double z_max = range[k][1];
        double t_elapsed = time_evaluation(z_min, z_max);
        printf("[z_min,z_max] = [ %7.1f, %7.1f ]",z_min,z_max);
        printf(" t_elapsed = %8.2f ns\n",1.E9*t_elapsed);
    }
}
