// SPDX-FileCopyrightText: 2025 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep/fluds/aahd_structs.h"
#include <cstdint>

namespace opensn
{
class DiscreteOrdinatesProblem;
class LBSGroupset;
class AngleSet;
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
            AngleSet& angle_set,
            AAHD_FLUDS& fluds,
            bool is_surface_source_active);

  // mesh and quadrature
  const char* mesh_data;
  const char* quad_data;
  // source moments and phi
  const double* src_moment;
  double* phi;
  // angle set
  const std::uint32_t* directions;
  std::uint32_t angleset_size;
  // group set
  std::uint32_t num_groups;
  std::uint32_t groupset_start;
  std::uint32_t groupset_size;
  // fluds
  AAHD_FLUDSPointerSet flud_data;
  const std::uint64_t* flud_index;
};

} // namespace opensn::gpu_kernel
