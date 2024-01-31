#pragma once

#include "mpicpp-lite/mpicpp-lite.h"
#include <vector>

typedef std::vector<double> VecDbl;

namespace mpi = mpicpp_lite;

namespace opensn
{

/**
 * Computes a global L2-norm
 */
double Vec2NormMPI(const VecDbl& x, const mpi::Communicator& comm);

} // namespace opensn
