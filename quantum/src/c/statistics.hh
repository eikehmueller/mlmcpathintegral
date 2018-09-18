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
    reset();
  }

  /** @brief Return label */
  std::string label() const {
    return obj_label;
  }

  /** @brief Reset all counters */
  void reset() {
    n_samples = 0;
    Q_k.clear();
    S_k.clear();
    S_k.resize(k_max,0.0);
    avg = 0.0;
  }

  /** @brief Record a new sample
   * 
   * @param[in] Q Value of new sample
   */
  void record_sample(const double Q) {
    n_samples++;
    Q_k.push_front(Q);
    // Push out last sample from deque
    if (Q_k.size()>k_max) {
      Q_k.pop_back();
    }
    // Update running average
    avg = ((n_samples-1.0)*avg + Q)/(1.0*n_samples);
    // Update running S_k
    for (int k=0;k<Q_k.size();++k) {
      unsigned int N_k = n_samples - k;
      S_k[k] = ((N_k-1.0)*S_k[k] + Q_k[0]*Q_k[k])/(1.0*N_k);
    }
  }
  

  /** @brief Return estimator for variance
   */
  double variance() const {
    return 1.0*n_samples/(n_samples-1.0)*(S_k[0]-avg*avg);
  }

  /** @brief Return estimator for average */
  double average() const {
    return avg;
  }
  /** @brief Return estimator for error of average */
  double error() const {
    return sqrt(tau_int()*variance()/(1.0*n_samples));
  }

  /** @brief Return vector with autocorrelation function \f$C(k)\f$ */
  std::vector<double> auto_corr() const {
    std::vector<double> S_k_tmp(S_k.size());
    std::copy(S_k.begin(),S_k.end(),S_k_tmp.begin());
    for (auto it=S_k_tmp.begin();it!=S_k_tmp.end();++it) {
      *it -= avg*avg;
    }
    return S_k_tmp;
  }

  /** @brief Return integrated autocorrelation time \f$\tau_{\text{int}}\f$ */
  double tau_int() const {
    double tau_int_tmp=0.0;
    for(int k=1;k<Q_k.size();++k) {
      tau_int_tmp += (1.-k/(1.0*n_samples))*(S_k[k]-avg*avg);
    }
    return 1.0+2.0*tau_int_tmp/(S_k[0]-avg*avg);
  }

  /** @brief Return the number of samples */
  unsigned int samples() const { return n_samples; }
  
private:
  /** @brief Label for identifying object */
  const std::string obj_label;
  /** @brief Window over which auto-correlations are measured */
  const unsigned int k_max;
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
  /** @brief Running average */
  double avg;
};

std::ostream& operator<<(std::ostream& os, const Statistics& stats);

#endif // STATISTICS_HH
