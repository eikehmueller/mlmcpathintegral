#include "mpi_wrapper.hh"

/* Initialise MPI */
void mpi_init() {
#ifdef USE_MPI
  if (not mpi_initialised) {
    MPI_Init(NULL, NULL);
    mpi_initialised=true;
  }
#endif // USE_MPI
}

/* Finalise MPI */
void mpi_finalize() {
#ifdef USE_MPI
  MPI_Finalize();
#endif // USE_MPI
}

/* Return MPI rank */
int mpi_comm_rank() {
  int comm_rank=0;
#ifdef USE_MPI
  MPI_Comm_rank(MPI_COMM_WORLD, &comm_rank);
#endif // USE_MPI
  return comm_rank;
}

/* Return total number of MPI ranks */
int mpi_comm_size() {
  int comm_size=1;
#ifdef USE_MPI
  MPI_Comm_size(MPI_COMM_WORLD, &comm_size);
#endif // USE_MPI
  return comm_size;
}

/* Return true on MPI master */
bool mpi_master() {
  return (mpi_comm_rank() == 0);
}

/* Parallel sum for scalar valued quantities */
double mpi_allreduce_sum(const double x) {
  double recv_data=x;
#ifdef USE_MPI
  MPI_Allreduce(&x, &recv_data, 1, MPI::DOUBLE, MPI_SUM, MPI_COMM_WORLD);
#endif // USE_MPI
  return recv_data;
}

/* Parallel sum for scalar valued quantities */
double mpi_allreduce_avg(const double x) {
  double recv_data=x;
#ifdef USE_MPI
  MPI_Allreduce(&x, &recv_data, 1, MPI::DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  recv_data /= mpi_comm_size();
#endif // USE_MPI
  return recv_data;
}

/* Parallel average for scalar valued quantities */
std::vector<double> mpi_allreduce_avg(const std::vector<double> x) {
  std::vector<double> result;
#ifdef USE_MPI
  int n_data = mpi_allreduce_min(x.size());
  double* send_data = (double*) malloc(n_data*sizeof(double));
  double* recv_data = (double*) malloc(n_data*sizeof(double));
  for (int j=0;j<n_data;++j) {
    send_data[j] = x[j];
  }
  MPI_Allreduce(send_data, recv_data, n_data,
                MPI::DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  double inv_nproc = 1.0/mpi_comm_size();
  for (int j=0;j<n_data;++j) {
    result.push_back(inv_nproc*recv_data[j]);
  }
  free(send_data);
  free(recv_data);
#else
  std::copy(x.begin(),x.end(),result.begin());
#endif // USE_MPI
  return result;
}

/* Parallel sum for scalar valued quantities */ 
int mpi_allreduce_sum(const int x) {
  int recv_data=x;
#ifdef USE_MPI
  MPI_Allreduce(&x, &recv_data, 1, MPI::INT, MPI_SUM, MPI_COMM_WORLD);
#endif // USE_MPI
  return recv_data;
}

/* Parallel sum for scalar valued quantities */ 
int mpi_allreduce_min(const int x) {
  int recv_data=x;
#ifdef USE_MPI
  MPI_Allreduce(&x, &recv_data, 1, MPI::INT, MPI_MIN, MPI_COMM_WORLD);
#endif // USE_MPI
  return recv_data;
}

/* Parallel sum for scalar valued quantities */ 
unsigned int mpi_allreduce_sum(const unsigned int x) {
  int recv_data=x;
#ifdef USE_MPI
  MPI_Allreduce(&x, &recv_data, 1, MPI::UNSIGNED, MPI_SUM, MPI_COMM_WORLD);
#endif // USE_MPI
  return recv_data;
}

/* Parallel logical and for scalar valued quantities */
bool mpi_allreduce_and(const bool x) {
  bool recv_data=x;
#ifdef USE_MPI
  MPI_Allreduce(&x, &recv_data, 1, MPI::BOOL, MPI_LAND, MPI_COMM_WORLD);
#endif // USE_MPI
  return recv_data;
}


/* Broadcast double value to all processes */
void mpi_bcast(double& x,const int root) {
#ifdef USE_MPI
  MPI_Bcast(&x,1,MPI::DOUBLE,root,MPI_COMM_WORLD);
#endif // USE_MPI
}

/* Broadcast integer value to all processes */
void mpi_bcast(int& x,const int root) {
#ifdef USE_MPI
  MPI_Bcast(&x,1,MPI::INT,root,MPI_COMM_WORLD);
#endif // USE_MPI
}

/* Broadcast bool value to all processes */
void mpi_bcast(bool& x,const int root) {
#ifdef USE_MPI
  MPI_Bcast(&x,1,MPI_LOGICAL,root,MPI_COMM_WORLD);
#endif // USE_MPI
}

/* Broadcast string value to all processes */
void mpi_bcast(std::string& x,const int root) {
#ifdef USE_MPI
  int n = x.length();
  mpi_bcast(n);
  char* char_array;
  char_array = (char*) malloc(sizeof(char)*n);
  if (mpi_master()) {
    for (int j=0;j<n;++j) {
      char_array[j] = x[j];
    }
  }
  MPI_Bcast(char_array,n,MPI::CHAR,root,MPI_COMM_WORLD);
  x.resize(n);
  for (int j=0;j<n;++j) {
    x[j] = char_array[j];
  }
  free(char_array);
#endif // USE_MPI
}

/* Scatter list of unsigned int values  to all processes */
void mpi_scatter(unsigned int* data, unsigned int& x) {
#ifdef USE_MPI
  int n_rank = mpi_comm_size();
  MPI_Scatter(data,1,MPI::UNSIGNED,
              &x,1,MPI::UNSIGNED,
              0,MPI_COMM_WORLD);
#else
  x = data[0];  
#endif // USE_MPI
}

/* Call mpi_finalize() and exit */
void mpi_exit(const int exit_code) {
  mpi_finalize();
  exit(exit_code);
}

/* Distribute number between processors */
unsigned int distribute_n(const unsigned int n) {
#ifdef USE_MPI
  unsigned int n_local;
  int n_proc = mpi_comm_size();
  unsigned int* n_local_list = (unsigned int*) malloc(n_proc);
  unsigned int n_floor = (unsigned int) floor(n/(1.0*n_proc));
  unsigned int n_overflow = n-n_floor*n_proc;
  for (int j=0;j<n_proc;++j) {
    n_local_list[j] = n_floor;
    if (j<n_overflow) n_local_list[j]++;
  }
  MPI_Scatter(n_local_list,1,MPI::UNSIGNED,
              &n_local,1,MPI::UNSIGNED,
              0,MPI_COMM_WORLD);
  free(n_local_list);
  return n_local;
#else
  return n;
#endif // USE_MPI
}
