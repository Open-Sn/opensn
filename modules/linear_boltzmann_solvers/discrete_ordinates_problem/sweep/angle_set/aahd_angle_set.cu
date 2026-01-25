// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep/angle_set/aahd_angle_set.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep_chunks/sweep_chunk.h"
#include "framework/logging/log.h"
#include "framework/runtime.h"
#include "caliper/cali.h"

namespace opensn
{

AAHD_AngleSet::AAHD_AngleSet(size_t id,
                             size_t num_groups,
                             const SPDS& spds,
                             std::shared_ptr<FLUDS>& fluds,
                             std::vector<size_t>& angle_indices,
                             std::map<uint64_t, std::shared_ptr<SweepBoundary>>& boundaries,
                             int maximum_message_size,
                             const MPICommunicatorSet& comm_set)
  : AngleSet(id, num_groups, spds, fluds, angle_indices, boundaries),
    async_comm_(*fluds, num_groups_, angle_indices.size(), maximum_message_size, comm_set),
    angle_indices_pinner_(angles_)
{
  stream_ = crb::Stream::create();
  angle_indices_pinner_.CopyToDevice();
}

void
AAHD_AngleSet::UpdateSweepDependencies(std::set<AngleSet*>& following_angle_sets)
{
  std::transform(following_angle_sets.begin(),
                 following_angle_sets.end(),
                 std::back_inserter(following_angle_sets_),
                 [](AngleSet* as) { return static_cast<AAHD_AngleSet*>(as); });
  for (auto* following_angle_set : following_angle_sets_)
  {
    ++(following_angle_set->num_dependencies_);
  }
}

void
AAHD_AngleSet::SetStartingLatch()
{
  starting_latch_ = std::make_unique<std::latch>(num_dependencies_);
}

void
AAHD_AngleSet::InitializeDelayedUpstreamData()
{
  async_comm_.InitializeDelayedUpstreamData();
}

AngleSetStatus
AAHD_AngleSet::AngleSetAdvance(SweepChunk& sweep_chunk, AngleSetStatus permission)
{
  CALI_CXX_MARK_SCOPE("AAHD_AngleSet::AngleSetAdvance");

  if (executed_)
    return AngleSetStatus::FINISHED;

  async_comm_.ReceiveUpstreamPsi(static_cast<int>(this->GetID()), stream_);
  starting_latch_->wait();
  async_comm_.InitializeLocalAndDownstreamBuffers(stream_);

  sweep_chunk.Sweep(*this);

  async_comm_.SendDownstreamPsi(static_cast<int>(this->GetID()), stream_);
  async_comm_.ClearLocalAndReceiveBuffers(stream_);
  async_comm_.ReceiveDelayedData(static_cast<int>(this->GetID()));

  for (auto& following_as : following_angle_sets_)
  {
    following_as->starting_latch_->count_down();
  }
  async_comm_.Wait();

  executed_ = true;
  return AngleSetStatus::FINISHED;
}

void
AAHD_AngleSet::ResetSweepBuffers()
{
  async_comm_.Reset();
  executed_ = false;
}

} // namespace opensn
