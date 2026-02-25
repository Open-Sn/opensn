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
                             const MPICommunicatorSet& comm_set)
  : CBC_AngleSet(id, num_groups, spds, fluds, angle_indices, boundaries, comm_set),
    stream_(crb::Stream::create()),
    device_angle_indices_(angles_.size(), stream_)
{
  crb::MemoryPinningManager angle_indices_pinner_(angles_);
  crb::copy(device_angle_indices_, angle_indices_pinner_, angles_.size(), 0, 0, stream_);
  // Set CBCD_FLUDS stream and asynchronously allocate storage for local psi
  auto* cbcd_fluds = std::static_pointer_cast<CBCD_FLUDS>(fluds_).get();
  cbcd_fluds->GetStream() = stream_;
  cbcd_fluds->AllocateLocalAndSavedPsi();
}

CBCD_AngleSet::~CBCD_AngleSet()
{
  device_angle_indices_.async_free(stream_);
}

} // namespace opensn