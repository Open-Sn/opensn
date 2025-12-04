// SPDX-FileCopyrightText: 2025 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep_chunks/cbc_gpu_kernel/arguments.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/discrete_ordinates_problem.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep/fluds/cbcd_fluds_common_data.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep/fluds/cbcd_fluds.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep/angle_set/angle_set.h"
#include "modules/linear_boltzmann_solvers/lbs_problem/groupset/lbs_groupset.h"
#include "modules/linear_boltzmann_solvers/lbs_problem/device/carrier/mesh_carrier.h"
#include "modules/linear_boltzmann_solvers/lbs_problem/device/carrier/quadrature_carrier.h"
#include "modules/linear_boltzmann_solvers/lbs_problem/device/memory_pinner.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep/sweep.h"

namespace opensn
{

cbc_gpu_kernel::Arguments::Arguments(DiscreteOrdinatesProblem& problem,
                                     const LBSGroupset& groupset,
                                     AngleSet& angle_set,
                                     CBCD_FLUDS& fluds,
                                     std::vector<Task*>& tasks)
{
  // Get mesh and quadrature data
  auto* mesh = reinterpret_cast<MeshCarrier*>(problem.GetCarrier(2));
  mesh_data = mesh->GetDevicePtr();
  auto* quadrature = reinterpret_cast<QuadratureCarrier*>(groupset.quad_carrier);
  quad_data = quadrature->GetDevicePtr();

  // Copy source moment and destination phi data to device
  auto* src = reinterpret_cast<MemoryPinner<double>*>(problem.GetPinner(0));
  src_moment = src->GetDevicePtr();
  auto* scalar_flux = reinterpret_cast<MemoryPinner<double>*>(problem.GetPinner(1));
  phi = scalar_flux->GetDevicePtr();

  // Copy angleset data to device
  auto* directions_num = reinterpret_cast<MemoryPinner<std::uint32_t>*>(angle_set.GetMemoryPin());
  directions = directions_num->GetDevicePtr();
  angleset_size = angle_set.GetNumAngles();

  // Copy groupset data to device
  groupset_size = groupset.groups.size();
  groupset_start = groupset.groups.front().id;
  num_groups = problem.GetGroups().size();

  // Copy ready cell local IDs to device, and retrieve CBCD_FLUDS pointers and node indices
  cell_local_ids = fluds.GetLocalCellIDs().GetDevicePtr();
  flud_data = fluds.GetPointerSet();
  flud_index = fluds.GetCommonData().GetDeviceCellFaceNodeMap();

  // Set batch size
  batch_size = static_cast<std::uint32_t>(tasks.size() * angleset_size * groupset_size);
}

} // namespace opensn