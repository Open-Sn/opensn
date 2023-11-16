#include "framework/math/math_mpi.h"

#include "framework/math/math.h"

namespace opensn
{

double
Vec2NormMPI(const VecDbl& x, MPI_Comm comm)
{
  size_t n = x.size();
  double local_sum = 0.;

  for (size_t i = 0; i != n; i++)
    local_sum += x[i] * x[i];

  double global_sum;
  MPI_Allreduce(&local_sum, &global_sum, 1, MPI_DOUBLE, MPI_SUM, comm);

  return sqrt(global_sum);
}

} // namespace opensn
