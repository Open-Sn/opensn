// SPDX-FileCopyrightText: 2025 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "modules/linear_boltzmann_solvers/lbs_problem/device/debug_tool/print_xs.h"

namespace opensn
{

__global__ void
PrintXS(double* xs, int block_id, std::uint32_t num_groups)
{
  std::printf("[DB]   Total XS for block ID %d at %p:", block_id, xs);
  for (std::uint32_t g = 0; g < num_groups; ++g)
  {
    std::printf(" %f", xs[g]);
  }
  std::printf("\n");
}

void
debug::PrintTotalXSOnDevice(TotalXSCarrier& xs)
{
  std::printf("[DB] Print total cross sections on GPU for each block ID:\n");
  for (auto& [block_id, index] : xs.block_id_to_index)
  {
    PrintXS<<<1, 1>>>(xs.GetXSGPUData(block_id), block_id, xs.num_groups);
    crb::synchronize();
  }
}

} // namespace opensn
