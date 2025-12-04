// SPDX-FileCopyrightText: 2026 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep/angle_set/cbc_angle_set.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep/fluds/cbcd_fluds.h"
#include "caribou/main.hpp"
#include <memory>

namespace crb = caribou;

namespace opensn
{

struct Task;

/**
 * \brief CBC angle set for device.
 */
class CBCD_AngleSet : public CBC_AngleSet
{
public:
  CBCD_AngleSet(size_t id,
                size_t num_groups,
                const SPDS& spds,
                std::shared_ptr<FLUDS>& fluds,
                const std::vector<size_t>& angle_indices,
                std::map<std::uint64_t, std::shared_ptr<SweepBoundary>>& boundaries,
                const MPICommunicatorSet& comm_set,
                bool use_gpu);

  /// Get the associated stream.
  crb::Stream& GetStream() { return stream_; }

  /// Get the current task list.
  std::vector<Task>& GetCurrentTaskList() { return current_task_list_; }

private:
  /// crb::Stream for this angle set.
  crb::Stream stream_;
};

} // namespace opensn