// SPDX-FileCopyrightText: 2025 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep_chunks/gpu_kernel/arguments.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/discrete_ordinates_problem.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep/fluds/aahd_fluds_common_data.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep/fluds/aahd_fluds.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep/fluds/cbcd_fluds_common_data.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep/fluds/cbcd_fluds.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep/angle_set/angle_set.h"
#include "modules/linear_boltzmann_solvers/lbs_problem/groupset/lbs_groupset.h"
#include "modules/linear_boltzmann_solvers/lbs_problem/device/carrier/mesh_carrier.h"
#include "modules/linear_boltzmann_solvers/lbs_problem/device/carrier/quadrature_carrier.h"
#include "modules/linear_boltzmann_solvers/lbs_problem/device/memory_pinner.h"

namespace opensn::gpu_kernel
{

Arguments::Arguments(DiscreteOrdinatesProblem& problem,
                     const LBSGroupset& groupset,
                     AngleSet& angle_set)
{
  // Get source moment and destination phi data
  auto* src = reinterpret_cast<MemoryPinner<double>*>(problem.GetPinner(0));
  src_moment = src->GetDevicePtr();
  MemoryPinner<double>* scalar_flux = reinterpret_cast<MemoryPinner<double>*>(problem.GetPinner(1));
  phi = scalar_flux->GetDevicePtr();
  // Get mesh and quadrature data
  auto* mesh = reinterpret_cast<MeshCarrier*>(problem.GetCarrier(2));
  mesh_data = mesh->GetDevicePtr();
  auto* quadrature = reinterpret_cast<QuadratureCarrier*>(groupset.quad_carrier);
  quad_data = quadrature->GetDevicePtr();
  // Copy angleset data to GPU
  auto* directions_num = reinterpret_cast<MemoryPinner<std::uint32_t>*>(angle_set.GetMemoryPin());
  directions = directions_num->GetDevicePtr();
  angleset_size = angle_set.GetNumAngles();
  // Copy groupset data to GPU
  groupset_size = groupset.groups.size();
  groupset_start = groupset.groups.front().id;
  num_groups = problem.GetGroups().size();
}

AAH_Arguments::AAH_Arguments(DiscreteOrdinatesProblem& problem,
                             const LBSGroupset& groupset,
                             AngleSet& angle_set,
                             AAHD_FLUDS& fluds,
                             bool is_surface_source_active)
  : Arguments(problem, groupset, angle_set)
{
  // Copy source moment and destination phi data to GPU
  auto* src = reinterpret_cast<MemoryPinner<double>*>(problem.GetPinner(0));
  src->CopyToDevice();
  MemoryPinner<double>* scalar_flux = reinterpret_cast<MemoryPinner<double>*>(problem.GetPinner(1));
  scalar_flux->CopyToDevice();
  // Copy FLUDS data to GPU and retrieve the pointer set
  flud_data =
    fluds.PrepareForSweep(*(problem.GetGrid()), angle_set, groupset, is_surface_source_active);
  flud_index = fluds.GetCommonData().GetDeviceIndex();
}

} // namespace opensn::gpu_kernel

namespace opensn::cbc_gpu_kernel
{

CBC_Arguments::CBC_Arguments(DiscreteOrdinatesProblem& problem,
                             const LBSGroupset& groupset,
                             AngleSet& angle_set,
                             CBCD_FLUDS& fluds,
                             const size_t num_ready_cells)
  : Arguments(problem, groupset, angle_set)
{
  // Get device pointer to ready cell local IDs
  cell_local_ids = fluds.GetLocalCellIDs().GetDevicePtr();
  // Retrieve CBCD_FLUDS pointers and node indices
  flud_data = fluds.GetPointerSet();
  flud_index = fluds.GetCommonData().GetDeviceIndex();
  // Set batch size
  batch_size = static_cast<std::uint32_t>(num_ready_cells * angleset_size * groupset_size);
}

} // namespace opensn::cbc_gpu_kernel
