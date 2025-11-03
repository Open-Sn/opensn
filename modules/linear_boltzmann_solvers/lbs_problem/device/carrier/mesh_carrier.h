// SPDX-FileCopyrightText: 2025 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include "modules/linear_boltzmann_solvers/lbs_problem/device/carrier/carrier.h"
#include "modules/linear_boltzmann_solvers/lbs_problem/device/carrier/outflow_carrier.h"
#include "modules/linear_boltzmann_solvers/lbs_problem/device/carrier/total_xs_carrier.h"
#include "modules/linear_boltzmann_solvers/lbs_problem/lbs_problem.h"
#include <cstdint>
#include <vector>

namespace opensn
{

/// Object managing mesh data on GPU (for transport only).
class MeshCarrier : public Carrier
{
public:
  /// Constructor.
  MeshCarrier(LBSProblem& lbs_problem, TotalXSCarrier& xs, OutflowCarrier& outflow);

  /// Total number of cell nodes in the mesh.
  std::uint64_t num_nodes_total;
  /// Vector storing offset for save angular flux for each cell.
  std::vector<std::uint64_t> saved_psi_offset;

protected:
  /// Compute the memory size (in bytes) to allocate on the GPU to copy the mesh data over.
  std::uint64_t ComputeSize(LBSProblem& lbs_problem);
  /**
   * Gather the data of each cell and faces into a contiguous chunk of memory on the host before
   * copying the data to the device.
   */
  void Assemble(LBSProblem& lbs_problem, TotalXSCarrier& xs, OutflowCarrier& outflow);
};

} // namespace opensn
