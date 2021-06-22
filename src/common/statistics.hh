#ifndef STATISTICS_HH
#define STATISTICS_HH STATISTICS_HH

/** @file statistics.hh
 * @brief Header file for statistics class
 */

#include <deque>
#include <vector>
#include <cmath>
#include <iostream>
#include <string>
#include <iomanip>
#include "common/parameters.hh"
#include "mpi/mpi_wrapper.hh"

/** @class StatisticsParameters
 *
 * @brief Class for storing parameters for recording statistics.
 */
class StatisticsParameters : public Parameters {
public:
    /** @brief Construct a new instance */
    StatisticsParameters() :
        Parameters("statistics"),
        n_autocorr_window_(20),
        n_min_samples_qoi_(100) {
        addKey("n_autocorr_window",Integer,Positive);
        addKey("n_min_samples_qoi",Integer,Positive);
    }

    /** @brief Read parameters from file
     *
     * @param[in] filename Name of file to read
     */
    int readFile(const std::string filename) {

        int readSuccess = Parameters::readFile(filename);
        if (!readSuccess) {
            n_autocorr_window_ = getContents("n_autocorr_window")->getInt();
            n_min_samples_qoi_ = getContents("n_min_samples_qoi")->getInt();
        }
        return readSuccess;
    }

    /** @brief Return size of autocorrelation window */
    unsigned int n_autocorr_window() const {
        return n_autocorr_window_;
    }
    /** @brief Return minimal number of samples for QoI */
    unsigned int n_min_samples_qoi() const {
        return n_min_samples_qoi_;
    }
private:
    /** @brief Size of autocorrelation window */
    unsigned int n_autocorr_window_;
    /** @brief Minimal number of samples for qoi */
    unsigned int n_min_samples_qoi_;
};


/** @class Statistics
 * @brief Class for recording statistics of an observable
 *
 * This class can be used to collect statistics on a random observable
 * \f$Q\f$. It allows calculation of the following quantities:
 *
 * - *Estimated average*
 *   \f[
 *     \overline{Q} = \langle Q_i \rangle = \frac{1}{N}\sum_{i=0}^{N-1} Q_i
 *   \f]
 *
 * - *Estimated variance*
 *   \f[
 *     \text{Var}[Q] = \frac{1}{N-1}\sum_{i=0}^{N-1}(Q_i-\overline{Q})^2
 *   \f]
 *
 * - *Autocorrelation function*
 *   \f[
 *     C(k) = \langle (Q_i-\overline{Q})(Q_{i-k}-\overline{Q}) \rangle\\
 *          = \langle Q_i Q_{i-k} \rangle - \overline{Q}^2
 *   \f]
 *
 * - *Integrated autocorrelation*
 *   \f[
 *     \tau_{\text{int}} = 1 + 2 \sum_{k=1}^{N-1}\left(1-\frac{k}{N}\right)
 *                         \frac{C(k)}{C(0)}
 *   \f]
 *
 * - *Error of average estimator*
 *   \f[
 *     \delta \overline{Q} = \sqrt{\tau_{\text{int}}\frac{\text{Var}[Q]}{N}}
 *   \f]
 *
 * To achieve this, the calculates the quantities
 * \f[
 *    S_k = \frac{1}{N_k} \sum_{i=k}^{N-1} S_i S_{i-k}
 * \f]
 * for \f$k=0,\dots,k_{\max}\f$ where \f$N_k:=N-k\f$.
 *
 * The calculation of autocorrelation time and variance will only be
 * reset if the hard_reset() method is called. This ensures that statistics
 * for those quantities can already be collected in the warmup-phase.
 *
 * The class is inherently parallel, i.e. it will always report the statistics
 * accumulated across all processors.
 *
 */
class Statistics {
public:
    /** @brief Create a new instance
     *
     * @param[in] label Label for identifying object
     * @param[in] k_max_ Window over which autocorrelations are measured.
     */
    Statistics(const std::string label_,
               const unsigned int k_max_) :
        obj_label(label_), k_max(k_max_) {
        hard_reset();
    }

    /** @brief Return label */
    std::string label() const {
        return obj_label;
    }

    /** @brief Reset all counters */
    void reset() {
        n_samples = 0;
        avg = 0.0;
    }

    /** @brief Reset all counters, including the ones for long term statistics */
    void hard_reset() {
        reset();
        Q_k.clear();
        S_k.clear();
        S_k.resize(k_max,0.0);
        avg_longterm = 0.0;
        avg2_longterm = 0.0;
        avg3_longterm = 0.0;
        avg4_longterm = 0.0;
        n_samples_longterm = 0;
    }


    /** @brief Record a new sample
     *
     * @param[in] Q Value of new sample
     */
    void record_sample(const double Q);

    /** @brief Return estimator for variance
     */
    double variance() const;
    
    /* @brief Return estimator for error of variance
     * 
     * The estimator for the MSE of the variance is given by
     * 
     * RMSE[Var(Q)] = sqrt( 1/n * ( E[(Q-Q[Q])^4] - (E[(Q-E[Q])^2])^2 ) )
     * 
     * This is a biased estimator, but the bias is expected to
     * be small as n >> 1.
     */
    double variance_error() const;
    
    /** @brief Return estimator for average */
    double average() const;

    /** @brief Return estimator for error of average */
    double error() const;

    /** @brief Return vector with autocorrelation function \f$C(k)\f$ */
    std::vector<double> auto_corr() const;

    /** @brief Return integrated autocorrelation time \f$\tau_{\text{int}}\f$ */
    double tau_int() const;
    
    /** @brief Return size of autocorrelation window */
    unsigned int autocorr_window() const {
        return k_max;
    }

    /** @brief Return the number of global samples */
    unsigned int samples() const;

    /** @brief Return number of local samples
     *
     * This returns the number of samples collected on each process
     */
    unsigned int local_samples() const;

private:
    /** @brief Label for identifying object */
    const std::string obj_label;
    /** @brief Window over which auto-correlations are measured */
    const unsigned int k_max;
    /** @brief Number of collected samples for long term averages */
    unsigned int n_samples_longterm;
    /** @brief Number of collected samples */
    unsigned int n_samples;
    /** @brief Deque holding the last samples
     * This is necessary to update autocorrelations. The deque stores (in this
     * order) \f$Q_j, Q_{j-1}, \dots, Q_{j-k_{\max}}\f$.
     */
    std::deque<double> Q_k;
    /** @brief Vector with estimated autocorrelations
     * This stores (in this order) \f$S_0,S_1,\dots,S_{k_{\max}}\f$.
     */
    std::vector<double> S_k;
    /** @brief Running average for quantity */
    double avg;
    /** @brief Running average for long term quantities */
    double avg_longterm;
    /** @brief Running average for Q^2 */
    double avg2_longterm;
    /** @brief Running average for Q^3 */
    double avg3_longterm;
    /** @brief Running average for Q^4 */
    double avg4_longterm;
};

std::ostream& operator<<(std::ostream& os, const Statistics& stats);

#endif // STATISTICS_HH
