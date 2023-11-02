#include "chi_math_mpi.h"

#include "chi_math.h"

namespace chi_math
{

/**Computes a global L2-norm*/
double Vec2NormMPI(const VecDbl& x, MPI_Comm comm)
{
  size_t n = x.size();
  double local_sum = 0.;

  for(size_t i = 0; i != n; i++)
    local_sum += x[i]*x[i];

  double global_sum;
  MPI_Allreduce(&local_sum,       //sendbuf
                &global_sum,      //recvbuf
                1, MPI_DOUBLE,     //count + datatype
                MPI_SUM,           //operation
                comm);             //communicator

  return sqrt(global_sum);
}

}//namespace chi_math