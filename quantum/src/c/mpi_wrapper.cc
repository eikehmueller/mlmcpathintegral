#include "mpi_wrapper.hh"

/* Initialise MPI */
void mpi_init() {
# ifdef USE_MPI
  MPI_Init(NULL, NULL);
#endif // USE_MPI
}

/* Finalise MPI */
void mpi_finalize() {
# ifdef USE_MPI
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
double mpi_allreduce_sum_scalar(const double x) {
  double recv_data=x;
#ifdef USE_MPI
  MPI_Allreduce(&x, &recv_data, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD);
#endif // USE_MPI
  return recv_data;
}

/* Parallel sum for scalar valued quantities */ 
int mpi_allreduce_sum_scalar(const int x) {
  int recv_data=x;
#ifdef USE_MPI
  MPI_Allreduce(&x, &recv_data, 1, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD);
#endif // USE_MPI
  return recv_data;
}

/* Broadcast double value to all processes */
void mpi_bcast(double& x,const int root) {
#ifdef USE_MPI
  MPI_Bcast(&x,1,MPI_DOUBLE_PRECISION,root,MPI_COMM_WORLD);
#endif // USE_MPI
}

/* Broadcast integer value to all processes */
void mpi_bcast(int& x,const int root) {
#ifdef USE_MPI
  MPI_Bcast(&x,1,MPI_INTEGER,root,MPI_COMM_WORLD);
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
  MPI_Bcast(char_array,n,MPI_CHARACTER,root,MPI_COMM_WORLD);
  x.resize(n);
  for (int j=0;j<n;++j) {
    x[j] = char_array[j];
  }
  free(char_array);
#endif // USE_MPI
}
