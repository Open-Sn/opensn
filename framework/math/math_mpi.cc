#include "framework/math/math_mpi.h"

#include "framework/math/math.h"

namespace opensn
{

double
Vec2NormMPI(const VecDbl& x, const mpi::Communicator& comm)
{
  size_t n = x.size();
  double sum = 0.;

  for (size_t i = 0; i != n; i++)
    sum += x[i] * x[i];

  comm.all_reduce(sum, mpi::op::sum<double>());

  return sqrt(sum);
}

} // namespace opensn
