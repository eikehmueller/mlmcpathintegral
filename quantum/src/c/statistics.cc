#include "statistics.hh"

/* Record a new sample */
void Statistics::record_sample(const double Q) {
  n_samples++;
  n_samples_longterm++;
  Q_k.push_front(Q);
  // Push out last sample from deque
  if (Q_k.size()>k_max) {
    Q_k.pop_back();
  }
  // Update running averages
  avg = ((n_samples-1.0)*avg + Q)/(1.0*n_samples);
  avg_longterm = ((n_samples_longterm-1.0)*avg_longterm + Q)/(1.0*n_samples_longterm);
  // Update running S_k
  for (unsigned int k=0;k<Q_k.size();++k) {
    unsigned int N_k = n_samples_longterm - k;
    S_k[k] = ((N_k-1.0)*S_k[k] + Q_k[0]*Q_k[k])/(1.0*N_k);
  }
}
  
/* Return estimator for variance */
double Statistics::variance() const {
  double avg_ = mpi_allreduce_avg(avg_longterm);
  double avg2_ = mpi_allreduce_avg(S_k[0]);
  unsigned int n_samples_ = mpi_allreduce_sum(n_samples_longterm);
  return 1.0*n_samples_/(n_samples_-1.0)*(avg2_-avg_*avg_);
}

/* Return estimator for average */
double Statistics::average() const {
  if (mpi_comm_size() == 1) {
    return avg;
  } else {
    return mpi_allreduce_avg(avg); 
  }
}

/* Return estimator for error of average */
double Statistics::error() const {
  return sqrt(tau_int()*variance()/(1.0*samples()));
}

/* Return vector with autocorrelation function \f$C(k)\f$ */
std::vector<double> Statistics::auto_corr() const {
  std::vector<double> S_k_result;
  double avg_ = mpi_allreduce_avg(avg_longterm); 
  if (mpi_comm_size() == 1) {
    S_k_result.resize(S_k.size());
    std::copy(S_k.begin(),S_k.end(),S_k_result.begin());
  } else {
    std::vector<double> S_k_global = mpi_allreduce_avg(S_k);
    S_k_result.resize(S_k_global.size());
    std::copy(S_k_global.begin(),S_k_global.end(),S_k_result.begin());
  }
  for (auto it=S_k_result.begin();it!=S_k_result.end();++it) {
    *it -= avg_*avg_;
  }
  return S_k_result;
}

/* Return integrated autocorrelation time \f$\tau_{\text{int}}\f$ */
double Statistics::tau_int() const {
  std::vector<double> C_k_ = auto_corr();
  unsigned n_samples_ = mpi_allreduce_sum(n_samples_longterm);
  double tau_int_tmp=0.0;
  for(unsigned int k=1;k<C_k_.size();++k) {
    tau_int_tmp += (1.-k/(1.0*n_samples_))*C_k_[k];
  }
  return fmax(1.0,1.0+2.0*tau_int_tmp/C_k_[0]);
}

/* Return the number of samples (across all processors) */
unsigned int Statistics::samples() const {
  return mpi_allreduce_sum(n_samples);
}

/* Return number of local samples */
unsigned int Statistics::local_samples() const {
  return n_samples;
}

/* Output statistics to stream object */
std::ostream& operator<<(std::ostream& os, const Statistics& stats) {
  os << " ";
  os << std::setprecision(6) << std::fixed;
  os << stats.label() << ": Avg +/- Err = " << stats.average();
  os << " +/- " << stats.error() << std::endl;
  os << " " << stats.label() << ": Var         = " << stats.variance() << std::endl;
  os << std::setprecision(3) << std::fixed;
  os << " " << stats.label() << ": tau_{int}   = " << stats.tau_int() << std::endl;
  os << " " << stats.label() << ": # samples   = " << stats.samples() << std::endl;
  return os;
}
