// SPDX-FileCopyrightText: 2025 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/discrete_ordinates_problem.h"
#include "modules/linear_boltzmann_solvers/lbs_problem/device/device_vector_mirror.h"
#include "modules/linear_boltzmann_solvers/lbs_problem/outflow/outflow_carrier.h"
#include "framework/utils/caliper_scopes.h"

namespace opensn
{

void
DiscreteOrdinatesProblem::CopyPhiAndSrcToDevice(std::size_t flux_block)
{
  if (!use_gpus_)
    return;
  CaliperRegionScope cali_gpu_transfer_scope("GPUTransfer", CaliperGPUTransferScopeDepth());
  auto* src = GetSourceMomentsPinner();
  src->CopyBlockToDevice(flux_block);
  DeviceVectorMirror<double>* phi = GetPhiPinner();
  phi->CopyBlockToDevice(flux_block);
  outflow_carrier_->CopyToDevice();
}

void
DiscreteOrdinatesProblem::CopyPhiAndOutflowBackToHost(std::size_t flux_block)
{
  if (!use_gpus_)
    return;
  CaliperRegionScope cali_gpu_transfer_scope("GPUTransfer", CaliperGPUTransferScopeDepth());
  auto* phi = GetPhiPinner();
  phi->CopyBlockFromDevice(flux_block);
  outflow_carrier_->CopyFromDevice();
}

} // namespace opensn
