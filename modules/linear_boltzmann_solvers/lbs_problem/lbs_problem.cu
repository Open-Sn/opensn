// SPDX-FileCopyrightText: 2025 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "modules/linear_boltzmann_solvers/lbs_problem/lbs_problem.h"
#include "modules/linear_boltzmann_solvers/lbs_problem/device/carrier/mesh_carrier.h"
#include "modules/linear_boltzmann_solvers/lbs_problem/device/carrier/total_xs_carrier.h"
#include "modules/linear_boltzmann_solvers/lbs_problem/device/device_vector_mirror.h"
#include "modules/linear_boltzmann_solvers/lbs_problem/outflow/outflow_carrier.h"
#include "framework/utils/error.h"

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
LBSProblem::CheckCapableDevices()
{
  std::uint32_t num_gpus = crb::get_num_gpus();
  OpenSnLogicalErrorIf(num_gpus == 0, "LBSProblem::CheckCapableDevices: No GPU detected.");
}

} // namespace opensn
