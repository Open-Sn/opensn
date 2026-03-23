// SPDX-FileCopyrightText: 2025 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep_chunks/gpu_kernel/arguments.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/discrete_ordinates_problem.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep/fluds/aahd_fluds_common_data.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep/fluds/aahd_fluds.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep/angle_set/aahd_angle_set.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep/fluds/cbcd_fluds_common_data.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep/fluds/cbcd_fluds.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep/angle_set/cbcd_angle_set.h"
#include "modules/linear_boltzmann_solvers/lbs_problem/groupset/lbs_groupset.h"
#include "modules/linear_boltzmann_solvers/lbs_problem/device/carrier/mesh_carrier.h"
#include "modules/linear_boltzmann_solvers/lbs_problem/device/carrier/quadrature_carrier.h"
#include "modules/linear_boltzmann_solvers/lbs_problem/device/memory_pinner.h"

namespace opensn::gpu_kernel
{

template <SweepType t>
Arguments<t>::Arguments(DiscreteOrdinatesProblem& problem,
                        const LBSGroupset& groupset,
                        AngleSetType& angle_set,
                        FLUDSType& fluds)
{
  // get mesh and quadrature data
  auto* mesh = problem.GetMeshCarrier();
  mesh_data = mesh->GetDevicePtr();
  auto* quadrature = reinterpret_cast<QuadratureCarrier*>(groupset.quad_carrier);
  quad_data = quadrature->GetDevicePtr();
  // copy source moment and destination phi data to GPU
  auto* src = problem.GetSourceMomentsPinner();
  src_moment = src->GetDevicePtr();
  MemoryPinner<double>* scalar_flux = problem.GetPhiPinner();
  phi = scalar_flux->GetDevicePtr();
  // Copy groupset data to GPU
  groupset_size = groupset.GetNumGroups();
  groupset_start = groupset.first_group;
  num_groups = problem.GetNumGroups();
  // Copy angleset data to GPU
  directions = angle_set.GetDeviceAngleIndices();
  angleset_size = angle_set.GetNumAngles();
  // Copy FLUDS data to GPU and retrieve the pointer set
  flud_data = fluds.GetDevicePointerSet();
  flud_index = fluds.GetCommonData().GetDeviceIndex();
}

template struct Arguments<SweepType::AAH>;
template struct Arguments<SweepType::CBC>;

} // namespace opensn::gpu_kernel
