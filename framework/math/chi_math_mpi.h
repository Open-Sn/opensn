#pragma once

#include "chi_mpi.h"
#include <vector>
typedef std::vector<double> VecDbl;

namespace chi_math
{
// 02 Vector operations
double Vec2NormMPI(const VecDbl& x, MPI_Comm comm);
} // namespace chi_math


