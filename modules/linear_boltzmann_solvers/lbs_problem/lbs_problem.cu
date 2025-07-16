// SPDX-FileCopyrightText: 2025 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "modules/linear_boltzmann_solvers/lbs_problem/lbs_problem.h"
#include "modules/linear_boltzmann_solvers/lbs_problem/device/memory_pinner.h"
#include "modules/linear_boltzmann_solvers/lbs_problem/device/carrier/mesh_carrier.h"
#include "modules/linear_boltzmann_solvers/lbs_problem/device/carrier/outflow_carrier.h"
#include "modules/linear_boltzmann_solvers/lbs_problem/device/carrier/total_xs_carrier.h"

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
  TotalXSCarrier* xs_ptr = new TotalXSCarrier(*this);
  OutflowCarrier* of_ptr = new OutflowCarrier(*this);
  MeshCarrier* msh_ptr = new MeshCarrier(*this, *xs_ptr, *of_ptr);
  carriers_[0] = xs_ptr;
  carriers_[1] = of_ptr;
  carriers_[2] = msh_ptr;
  // initialize pinners
  MemoryPinner<double>* src = new MemoryPinner<double>(q_moments_local_);
  MemoryPinner<double>* phi = new MemoryPinner<double>(phi_new_local_);
  pinners_[0] = src;
  pinners_[1] = phi;
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
  if (carriers_[0])
  {
    TotalXSCarrier* xs_ptr = reinterpret_cast<TotalXSCarrier*>(carriers_[0]);
    delete xs_ptr;
    carriers_[0] = nullptr;
  }
  if (carriers_[1])
  {
    OutflowCarrier* of_ptr = reinterpret_cast<OutflowCarrier*>(carriers_[1]);
    delete of_ptr;
    carriers_[1] = nullptr;
  }
  if (carriers_[2])
  {
    MeshCarrier* msh_ptr = reinterpret_cast<MeshCarrier*>(carriers_[2]);
    delete msh_ptr;
    carriers_[2] = nullptr;
  }
  // delete pinners
  if (pinners_[0])
  {
    MemoryPinner<double>* src = reinterpret_cast<MemoryPinner<double>*>(pinners_[0]);
    delete src;
    pinners_[0] = nullptr;
  }
  if (pinners_[1])
  {
    MemoryPinner<double>* phi = reinterpret_cast<MemoryPinner<double>*>(pinners_[1]);
    delete phi;
    pinners_[1] = nullptr;
  }
}

void
LBSProblem::CheckCapableDevices()
{
  std::uint32_t num_gpus = crb::get_num_gpus();
  if (num_gpus == 0)
  {
    throw std::runtime_error("No GPU detected.\n");
  }
}

} // namespace opensn
