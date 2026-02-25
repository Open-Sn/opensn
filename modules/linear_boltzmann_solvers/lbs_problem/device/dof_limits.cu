// SPDX-FileCopyrightText: 2026 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "modules/linear_boltzmann_solvers/lbs_problem/device/dof_limits.h"
#include "caribou/main.hpp"

namespace crb = caribou;

namespace opensn
{

std::uint32_t
ComputeMaxDofSharedMem()
{
  std::uint32_t shared_mem_size = crb::get_max_shared_memory_per_block();
  std::uint32_t max_dof = max_dof_gpu_register;
  std::uint32_t required_shared_mem = crb::num_cores_per_sm * max_dof * max_dof * sizeof(double);
  while (required_shared_mem <= shared_mem_size)
  {
    ++max_dof;
    required_shared_mem = crb::num_cores_per_sm * max_dof * max_dof * sizeof(double);
  }
  return max_dof;
}

std::uint32_t max_dof_gpu_shared_mem = max_dof_gpu_register;

} // namespace opensn
