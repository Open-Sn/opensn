// SPDX-FileCopyrightText: 2025 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include "modules/linear_boltzmann_solvers/lbs_problem/device/carrier/carrier.h"
#include "modules/linear_boltzmann_solvers/lbs_problem/lbs_problem.h"
#include <cstdint>
#include <limits>
#include <map>

namespace opensn
{

/**
 * \brief Object managing the total cross section map on GPU.
 * \details For each block ID, group-wise total cross sections are stored side-by-side in device
 * memory. The constructor initializes a map to record the stride associated with each block ID.
 */
class TotalXSCarrier : public Carrier
{
public:
  /// Constructor from a fully initialized LBS problem.
  TotalXSCarrier(LBSProblem& lbs_problem);

  /// Retrive the pointer to the total cross section data on GPU for a given block ID.
  double* GetXSGPUData(int block_id);

  /// Number of groups.
  std::size_t num_groups = std::numeric_limits<std::size_t>::max();
  /// Number of block IDs.
  std::uint32_t num_block_ids;
  /// Map from block ID to flatten index.
  std::map<int, std::uint32_t> block_id_to_index;

protected:
  /// Compute the memory size (in bytes) to allocate on the GPU to copy the cross section data over.
  std::uint64_t ComputeSize(LBSProblem& lbs_problem);
  /**
   * Gather the total cross sections into a contiguous chunk of memory on the host before copying
   * the data to the device.
   */
  void Assemble(LBSProblem& lbs_problem);
};

} // namespace opensn
