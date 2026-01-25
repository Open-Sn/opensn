// SPDX-FileCopyrightText: 2025 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep_chunks/gpu_kernel/arguments.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/discrete_ordinates_problem.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep/fluds/aahd_fluds_common_data.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep/fluds/aahd_fluds.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep/angle_set/angle_set.h"
#include "modules/linear_boltzmann_solvers/lbs_problem/groupset/lbs_groupset.h"
#include "modules/linear_boltzmann_solvers/lbs_problem/device/carrier/mesh_carrier.h"
#include "modules/linear_boltzmann_solvers/lbs_problem/device/carrier/quadrature_carrier.h"
#include "modules/linear_boltzmann_solvers/lbs_problem/device/memory_pinner.h"

namespace opensn
{

gpu_kernel::Arguments::Arguments(DiscreteOrdinatesProblem& problem,
                                 const LBSGroupset& groupset,
                                 AngleSet& angle_set,
                                 AAHD_FLUDS& fluds,
                                 bool is_surface_source_active)
{
  // get mesh and quadrature data
  auto* mesh = reinterpret_cast<MeshCarrier*>(problem.GetCarrier(2));
  mesh_data = mesh->GetDevicePtr();
  auto* quadrature = reinterpret_cast<QuadratureCarrier*>(groupset.quad_carrier);
  quad_data = quadrature->GetDevicePtr();
  // copy source moment and destination phi data to GPU
  auto* src = reinterpret_cast<MemoryPinner<double>*>(problem.GetPinner(0));
  src_moment = src->GetDevicePtr();
  MemoryPinner<double>* scalar_flux = reinterpret_cast<MemoryPinner<double>*>(problem.GetPinner(1));
  phi = scalar_flux->GetDevicePtr();
  // copy angleset data to GPU
  auto* directions_num = reinterpret_cast<MemoryPinner<std::uint32_t>*>(angle_set.GetMemoryPin());
  directions = directions_num->GetDevicePtr();
  angleset_size = angle_set.GetNumAngles();
  // copy groupset data to GPU
  groupset_size = groupset.groups.size();
  groupset_start = groupset.groups.front().id;
  num_groups = problem.GetGroups().size();
  // copy FLUDS data to GPU and retrieve the pointer set
  flud_data =
    fluds.PrepareForSweep(*(problem.GetGrid()), angle_set, groupset, is_surface_source_active);
  flud_index = fluds.GetCommonData().GetDeviceIndex();
}

} // namespace opensn
