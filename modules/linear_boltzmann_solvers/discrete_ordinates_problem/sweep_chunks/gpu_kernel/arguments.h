// SPDX-FileCopyrightText: 2025 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep/fluds/aahd_structs.h"
#include <cstdint>

namespace opensn
{
class DiscreteOrdinatesProblem;
class LBSGroupset;
class AAHD_AngleSet;
class AAHD_FLUDS;
} // namespace opensn

namespace opensn::gpu_kernel
{

/// Arguments for the sweep kernel.
struct Arguments
{
  /// Constructor
  Arguments(DiscreteOrdinatesProblem& problem,
            const LBSGroupset& groupset,
            AAHD_AngleSet& angle_set,
            AAHD_FLUDS& fluds,
            bool is_surface_source_active);

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
  // fluds
  AAHD_FLUDSPointerSet flud_data;
  const std::uint64_t* __restrict__ flud_index;
};

} // namespace opensn::gpu_kernel
