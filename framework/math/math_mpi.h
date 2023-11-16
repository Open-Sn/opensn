#pragma once

#include "framework/mpi/mpi.h"
#include <vector>
typedef std::vector<double> VecDbl;

namespace opensn
{

/**
 * Computes a global L2-norm
 */
double Vec2NormMPI(const VecDbl& x, MPI_Comm comm);

} // namespace opensn
