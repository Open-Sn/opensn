// SPDX-FileCopyrightText: 2026 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep/angle_set/cbcd_angle_set.h"

namespace opensn
{

CBCD_AngleSet::CBCD_AngleSet(size_t id,
                             size_t num_groups,
                             const SPDS& spds,
                             std::shared_ptr<FLUDS>& fluds,
                             const std::vector<size_t>& angle_indices,
                             std::map<uint64_t, std::shared_ptr<SweepBoundary>>& boundaries,
                             const MPICommunicatorSet& comm_set,
                             bool use_gpu)
  : CBC_AngleSet(id, num_groups, spds, fluds, angle_indices, boundaries, comm_set, use_gpu),
    stream_(crb::Stream::create())
{
  std::static_pointer_cast<CBCD_FLUDS>(fluds_)->GetStream() = stream_;
}

} // namespace opensn