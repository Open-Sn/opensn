// SPDX-FileCopyrightText: 2025 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "modules/linear_boltzmann_solvers/lbs_problem/lbs_problem.h"
#include "modules/linear_boltzmann_solvers/lbs_problem/device/carrier/mesh_carrier.h"
#include "modules/linear_boltzmann_solvers/lbs_problem/device/carrier/total_xs_carrier.h"
#include "modules/linear_boltzmann_solvers/lbs_problem/device/device_vector_mirror.h"
#include "modules/linear_boltzmann_solvers/lbs_problem/outflow/outflow_carrier.h"
#include "framework/utils/error.h"
#include "framework/utils/caliper_scopes.h"

namespace opensn
{

void
LBSProblem::InitializeGPUExtras()
{
  // exit if GPU acceleration is not enabled
  if (!use_gpus_)
  {
    return;
  }
  CaliperRegionScope cali_gpu_setup_scope("GPUSetup", CaliperGPUSetupScopeDepth());
  // initialize carriers
  total_xs_carrier_ = std::make_shared<TotalXSCarrier>(*this);
  outflow_carrier_ = std::make_shared<OutflowCarrier>(*this);
  mesh_carrier_ = std::make_shared<MeshCarrier>(*this, *total_xs_carrier_, *outflow_carrier_);
  // initialize pinners
  source_pinner_ = std::make_shared<DeviceVectorMirror<double>>(q_moments_local_);
  phi_pinner_ = std::make_shared<DeviceVectorMirror<double>>(phi_new_local_);
}

void
LBSProblem::ResetGPUCarriers()
{
  // exit if GPU acceleration is not enabled
  if (!use_gpus_)
  {
    return;
  }
  // delete carriers
  total_xs_carrier_.reset();
  outflow_carrier_.reset();
  mesh_carrier_.reset();
  // delete pinners
  source_pinner_.reset();
  phi_pinner_.reset();
}

void
LBSProblem::RebuildOutflowDependentGPUCarriers()
{
  // exit if GPU acceleration is not enabled
  if (not use_gpus_)
    return;
  // outflow_carrier_ gets a fresh device allocation (a different address) every time it's
  // rebuilt, and mesh_carrier_'s own GPU buffer bakes in raw outflow_carrier_ device pointers
  // (see MeshCarrier::Assemble) -- both must be rebuilt together, or mesh_carrier_ is left
  // pointing at freed device memory. total_xs_carrier_'s device pointers are ALSO baked into
  // mesh_carrier_ the same way, but total_xs_carrier_ itself is untouched here, so its address
  // (and therefore what mesh_carrier_ re-reads from it) doesn't change -- safe to leave alone,
  // along with source_pinner_/phi_pinner_, which have no outflow/XS dependency at all.
  //
  // Unlike ResetGPUCarriers+InitializeGPUExtras, this relies on total_xs_carrier_ already being
  // built (it's read, not rebuilt, above) -- fall back to a full rebuild rather than
  // dereferencing a null carrier if that precondition somehow doesn't hold, e.g. if this is ever
  // called before the normal InitializeGPUExtras() call in LBSProblem::InitializeRuntimeCore has
  // run.
  if (not total_xs_carrier_)
  {
    ResetGPUCarriers();
    InitializeGPUExtras();
    return;
  }
  outflow_carrier_.reset();
  mesh_carrier_.reset();
  outflow_carrier_ = std::make_shared<OutflowCarrier>(*this);
  mesh_carrier_ = std::make_shared<MeshCarrier>(*this, *total_xs_carrier_, *outflow_carrier_);
}

void
LBSProblem::CheckCapableDevices()
{
  std::uint32_t num_gpus = crb::get_num_gpus();
  OpenSnLogicalErrorIf(num_gpus == 0, "LBSProblem::CheckCapableDevices: No GPU detected.");
}

} // namespace opensn
