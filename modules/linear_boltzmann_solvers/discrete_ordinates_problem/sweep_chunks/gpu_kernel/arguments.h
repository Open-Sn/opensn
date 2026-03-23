// SPDX-FileCopyrightText: 2025 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep/fluds/aahd_structs.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep/fluds/cbcd_structs.h"
#include "caribou/main.hpp"
#include <cstdint>
#include <type_traits>

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

enum class SweepType
{
  AAH = 0,
  CBC = 1
};

consteval bool
to_bool(SweepType t)
{
  return t == SweepType::AAH;
}

template <SweepType t>
using NodeIndexType = std::conditional_t<to_bool(t), AAHD_NodeIndex, CBCD_NodeIndex>;

/// Argument structure for both AAH and CBC kernels
template <SweepType t>
struct Arguments
{
  using AngleSetType = std::conditional_t<to_bool(t), AAHD_AngleSet, CBCD_AngleSet>;
  using FLUDSType = std::conditional_t<to_bool(t), AAHD_FLUDS, CBCD_FLUDS>;
  using FLUDSPointerSetType =
    std::conditional_t<to_bool(t), AAHD_FLUDSPointerSet, CBCD_FLUDSPointerSet>;

  Arguments(DiscreteOrdinatesProblem& problem,
            const LBSGroupset& groupset,
            AngleSetType& angle_set,
            FLUDSType& fluds);

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
  FLUDSPointerSetType flud_data;
};

} // namespace opensn::gpu_kernel
