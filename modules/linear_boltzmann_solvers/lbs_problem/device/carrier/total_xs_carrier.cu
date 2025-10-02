// SPDX-FileCopyrightText: 2025 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "modules/linear_boltzmann_solvers/lbs_problem/device/carrier/total_xs_carrier.h"

namespace opensn
{

TotalXSCarrier::TotalXSCarrier(LBSProblem& lbs_problem)
{
  std::uint64_t size = ComputeSize(lbs_problem);
  host_memory_.reserve(size);
  host_memory_.resize(size);
  Assemble(lbs_problem);
  device_memory_ = crb::DeviceMemory<char>(size);
  crb::copy(device_memory_, host_memory_, size);
}

std::uint64_t
TotalXSCarrier::ComputeSize(LBSProblem& lbs_problem)
{
  std::uint64_t alloc_size = 0;
  const std::map<int, std::shared_ptr<MultiGroupXS>>& xs_map = lbs_problem.GetMatID2XSMap();
  // check if all cross sections have the same number of group
  for (const auto& [block_id, xs] : xs_map)
  {
    if (num_groups == std::numeric_limits<size_t>::max())
    {
      num_groups = xs->GetNumGroups();
    }
    else if (num_groups != xs->GetNumGroups())
    {
      throw std::runtime_error("Provided cross sections don't have the same number of groups.\n");
    }
  }
  // compute size
  num_block_ids = xs_map.size();
  alloc_size += num_block_ids * num_groups * sizeof(double);
  return alloc_size;
}

void
TotalXSCarrier::Assemble(LBSProblem& lbs_problem)
{
  char* data = reinterpret_cast<char*>(host_memory_.data());
  const std::map<int, std::shared_ptr<MultiGroupXS>>& xs_map = lbs_problem.GetMatID2XSMap();
  // copy total cross section data
  for (std::uint64_t index = 0; const auto& [block_id, xs] : xs_map)
  {
    block_id_to_index[block_id] = index++;
    double* total_xs_data = reinterpret_cast<double*>(data);
    const std::vector<double>& sigma_total = xs->GetSigmaTotal();
    std::copy(sigma_total.begin(), sigma_total.end(), total_xs_data);
    data = reinterpret_cast<char*>(total_xs_data + sigma_total.size());
  }
}

double*
TotalXSCarrier::GetXSGPUData(int block_id)
{
  char* gpu_data = device_memory_.get();
  return reinterpret_cast<double*>(gpu_data) + block_id_to_index[block_id] * num_groups;
}

} // namespace opensn
