// SPDX-FileCopyrightText: 2025 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/discrete_ordinates_problem.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep/fluds/aahd_structs.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep/fluds/aahd_fluds_common_data.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep/fluds/aahd_fluds.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep/angle_set/aahd_angle_set.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep/fluds/cbcd_fluds_common_data.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep/fluds/cbcd_fluds.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep/fluds/cbcd_structs.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep/angle_set/cbcd_angle_set.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep/boundary/boundary_carrier.h"
#include "modules/linear_boltzmann_solvers/lbs_problem/groupset/lbs_groupset.h"
#include "modules/linear_boltzmann_solvers/lbs_problem/device/carrier/mesh_carrier.h"
#include "modules/linear_boltzmann_solvers/lbs_problem/device/carrier/quadrature_carrier.h"
#include "modules/linear_boltzmann_solvers/lbs_problem/device/device_vector_mirror.h"
#include "caribou/main.hpp"
#include <cstdint>
#include <type_traits>

namespace crb = caribou;

namespace opensn::gpu_kernel
{

consteval bool
to_bool(SweepKind k)
{
  return k == SweepKind::AAH;
}

template <SweepKind k>
using NodeIndexType = std::conditional_t<to_bool(k), AAHD_NodeIndex, CBCD_NodeIndex>;

/// Arguments for AAHD and CBCD kernels
template <SweepKind k>
struct Arguments
{
  using AngleSetType = std::conditional_t<to_bool(k), AAHD_AngleSet, CBCD_AngleSet>;
  using FLUDSType = std::conditional_t<to_bool(k), AAHD_FLUDS, CBCD_FLUDS>;
  using FLUDSPointerSetType =
    std::conditional_t<to_bool(k), AAHD_FLUDSPointerSet, CBCD_FLUDSPointerSet>;

  Arguments(DiscreteOrdinatesProblem& problem,
            const LBSGroupset& groupset,
            AngleSetType& angle_set,
            FLUDSType& fluds,
            bool is_active)
    : boundary(problem.GetBoundaryCarrier()->GetDevicePtr(groupset.id)),
      directions(angle_set.GetDeviceAngleIndices()),
      angleset_size(angle_set.GetNumAngles()),
      num_groups(problem.GetNumGroups()),
      groupset_start(groupset.first_group),
      groupset_size(groupset.GetNumGroups()),
      is_surface_source_active(is_active)
  {
    // Get mesh and quadrature data
    auto* mesh = problem.GetMeshCarrier();
    mesh_data = mesh->GetDevicePtr();
    auto quadrature = groupset.quad_carrier;
    quad_data = quadrature->GetDevicePtr();
    // Copy source moment and destination phi data to GPU
    auto* src = problem.GetSourceMomentsPinner();
    src_moment = src->GetDevicePtr();
    DeviceVectorMirror<double>* scalar_flux = problem.GetPhiPinner();
    phi = scalar_flux->GetDevicePtr();
    // Copy groupset data to GPU
    if constexpr (k == SweepKind::AAH)
      boundary_offset = angle_set.GetDeviceBoudnaryOffset();
    // Copy FLUDS data to GPU and retrieve the pointer set
    flud_data = fluds.GetDevicePointerSet();
    flud_index = fluds.GetCommonData().GetDeviceIndex();
    // Copy surface source active
  }

  // Mesh and quadrature
  const char* __restrict__ mesh_data = nullptr;
  const char* __restrict__ quad_data = nullptr;
  // Source moments and phi
  const double* __restrict__ src_moment = nullptr;
  double* __restrict__ phi = nullptr;
  // Boundary
  double* __restrict__ boundary = nullptr;
  // Angle set
  const std::uint32_t* __restrict__ directions = nullptr;
  const std::uint64_t* __restrict__ boundary_offset = nullptr;
  std::uint32_t angleset_size = 0;
  // Group set
  std::uint32_t num_groups = 0;
  std::uint32_t groupset_start = 0;
  std::uint32_t groupset_size = 0;
  // FLUDS
  const std::uint64_t* __restrict__ flud_index = nullptr;
  FLUDSPointerSetType flud_data;
  // Source active
  bool is_surface_source_active = false;
};

} // namespace opensn::gpu_kernel
