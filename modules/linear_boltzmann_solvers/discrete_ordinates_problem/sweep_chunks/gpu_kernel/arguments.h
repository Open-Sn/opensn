// SPDX-FileCopyrightText: 2025 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep/fluds/aahd_structs.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep/fluds/cbcd_structs.h"
#include <cstdint>

namespace opensn
{
class DiscreteOrdinatesProblem;
class LBSGroupset;
class AngleSet;
class AAHD_FLUDS;
class CBCD_FLUDS;
} // namespace opensn

namespace opensn::gpu_kernel
{

/// Base arguments containing common data for both AAH and CBC kernels.
struct Arguments
{
  Arguments(DiscreteOrdinatesProblem& problem, const LBSGroupset& groupset, AngleSet& angle_set);

  // mesh and quadrature
  const char* __restrict__ mesh_data;
  const char* __restrict__ quad_data;
  // source moments and phi
  const double* __restrict__ src_moment;
  double* __restrict__ phi;
  // angle set
  const std::uint32_t* __restrict__ directions;
  std::uint32_t angleset_size;
  // group set
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
                AngleSet& angle_set,
                AAHD_FLUDS& fluds,
                bool is_surface_source_active);

  // FLUDS
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
                AngleSet& angle_set,
                CBCD_FLUDS& fluds,
                const size_t num_ready_cells);

  // Cell local IDs corresponding to the current set of ready cells that can be swept
  const std::uint64_t* cell_local_ids;
  // Device CBC_FLUDS pointers
  CBCD_FLUDSPointerSet flud_data;
  // Batch size info
  std::uint32_t batch_size;
};

} // namespace opensn::cbc_gpu_kernel
