// SPDX-FileCopyrightText: 2025 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include "modules/linear_boltzmann_solvers/lbs_problem/device/carrier/carrier.h"
#include "modules/linear_boltzmann_solvers/lbs_problem/lbs_problem.h"
#include <cstdint>
#include <map>

namespace opensn
{

/**
 * @brief Object managing outflow data on GPU.
 * @details This class manages memory for storing group-wise outflow data for each boundary face and
 * provides cleaner way to transfer the data back to the host.
 */
class OutflowCarrier : public Carrier
{
public:
  /// Constructor.
  OutflowCarrier(LBSProblem& lbs_problem);

  /// Accumulate computed outflow on GPU back to CPU.
  void AccumulateBack(std::vector<CellLBSView>& cell_transport_views);
  /**
   * Zero-fill device and host memory, forcing flux accumulation to zero to reuse the allocated
   * space for the next sweep.
   */
  void Reset(void);
  /// Get offset from cell local index and face index.
  std::uint64_t GetOffset(const std::uint32_t& cell_local_idx, const std::uint32_t& face_idx);

  /**
   * @brief Maps boundary faces to their outflow memory offsets.
   * @details This map associates each boundary face, which is identified by a tuple of cell local
   * index and face index in the cell, with its corresponding offset in the contiguous outflow data
   * stored on the device.
   *
   * @note Non-boundary faces are not included in this map.
   */
  std::map<std::uint64_t, std::uint64_t> outflow_map;
  /// Number of groups.
  std::uint64_t num_groups;
};

} // namespace opensn
