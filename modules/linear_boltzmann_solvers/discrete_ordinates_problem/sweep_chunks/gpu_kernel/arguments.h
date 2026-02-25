// SPDX-FileCopyrightText: 2025 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep/fluds/aahd_structs.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep/fluds/cbcd_structs.h"
#include "caribou/main.hpp"
#include <cstdint>

namespace crb = caribou;

namespace opensn
{
class DiscreteOrdinatesProblem;
class LBSGroupset;
class AAHD_AngleSet;
class AAHD_FLUDS;
class CBCD_AngleSet;
class CBCD_FLUDS;
} // namespace opensn

namespace opensn::gpu_kernel
{

#if defined(__NVCC__)
constexpr unsigned int threshold = 128;
#elif defined(__HIPCC__)
constexpr unsigned int threshold = 64;
#endif

static unsigned int
RoundUp(unsigned int num, unsigned int divisor = crb::get_warp_size())
{
  return (num + divisor - 1) & ~(divisor - 1);
}

/// Common arguments for AAH and CBC kernels.
struct Arguments
{
  Arguments(DiscreteOrdinatesProblem& problem, const LBSGroupset& groupset);

  // Mesh and quadrature
  const char* __restrict__ mesh_data;
  const char* __restrict__ quad_data;
  // Source moments and phi
  const double* __restrict__ src_moment;
  double* __restrict__ phi;
  // Angle set
  const std::uint32_t* __restrict__ directions;
  std::uint32_t angleset_size;
  // Group set
  std::uint32_t num_groups;
  std::uint32_t groupset_start;
  std::uint32_t groupset_size;
  // FLUDS
  const std::uint64_t* __restrict__ flud_index;
};

/// Arguments for the AAH sweep kernel.
struct AAH_Arguments : public Arguments
{
  AAH_Arguments(DiscreteOrdinatesProblem& problem,
                const LBSGroupset& groupset,
                AAHD_AngleSet& angle_set,
                AAHD_FLUDS& fluds,
                bool is_surface_source_active);

  // FLUDS pointer set
  AAHD_FLUDSPointerSet flud_data;
};

} // namespace opensn::gpu_kernel

namespace opensn::cbc_gpu_kernel
{

/// Arguments for the CBC sweep kernel.
struct CBC_Arguments : public opensn::gpu_kernel::Arguments
{
  CBC_Arguments(DiscreteOrdinatesProblem& problem,
                const LBSGroupset& groupset,
                CBCD_AngleSet& angle_set,
                CBCD_FLUDS& fluds);

  // Ready cell local IDs
  const std::uint64_t* cell_local_ids;
  // FLUDS pointer set
  CBCD_FLUDSPointerSet flud_data;
};

} // namespace opensn::cbc_gpu_kernel
