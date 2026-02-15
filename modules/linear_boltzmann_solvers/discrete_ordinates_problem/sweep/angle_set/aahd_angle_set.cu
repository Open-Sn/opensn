// SPDX-FileCopyrightText: 2026 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep/angle_set/aahd_angle_set.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep_chunks/aahd_sweep_chunk.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep/fluds/aahd_fluds.h"
#include "caliper/cali.h"
#include <algorithm>

namespace opensn
{

AAHD_AngleSet::AAHD_AngleSet(size_t id,
                             unsigned int num_groups,
                             const SPDS& spds,
                             std::shared_ptr<FLUDS>& fluds,
                             std::vector<size_t>& angle_indices,
                             std::map<uint64_t, std::shared_ptr<SweepBoundary>>& boundaries,
                             int maximum_message_size,
                             const MPICommunicatorSet& comm_set)
  : AngleSet(id, num_groups, spds, fluds, angle_indices, boundaries),
    async_comm_(*fluds, num_groups_, angle_indices.size(), maximum_message_size, comm_set),
    device_angle_indices_(angles_.size())
{
  stream_ = crb::Stream::create();
  std::dynamic_pointer_cast<AAHD_FLUDS>(fluds_)->GetStream() = stream_;
  crb::MemoryPinningManager angle_indices_pinner_(angles_);
  crb::copy(device_angle_indices_, angle_indices_pinner_, angles_.size());
  has_reflecting_boundaries_ =
    std::any_of(boundaries_.begin(),
                boundaries_.end(),
                [](auto& bndry_pair) { return bndry_pair.second->IsReflecting(); });
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
  auto* aahd_fluds = static_cast<AAHD_FLUDS*>(fluds_.get());
  aahd_fluds->AllocateDelayedPrelocIOutgoingPsi();
  aahd_fluds->AllocateDelayedLocalPsi();
}

AngleSetStatus
AAHD_AngleSet::AngleSetAdvance(SweepChunk& sweep_chunk, AngleSetStatus permission)
{
  CALI_CXX_MARK_SCOPE("AAHD_AngleSet::AngleSetAdvance");

  if (executed_)
    return AngleSetStatus::FINISHED;

  auto* aahd_fluds = static_cast<AAHD_FLUDS*>(fluds_.get());
  auto& aahd_sweep_chunk = static_cast<AAHDSweepChunk&>(sweep_chunk);

  async_comm_.PrepostReceiveUpstreamPsi(static_cast<int>(this->GetID()));
  async_comm_.PrepostReceiveDelayedData(static_cast<int>(this->GetID()));
  aahd_fluds->CopyDelayedPsiToDevice();
  starting_latch_->wait();

  aahd_fluds->CopyBoundaryToDevice(aahd_sweep_chunk.GetGrid(),
                                   *this,
                                   aahd_sweep_chunk.GetGroupset(),
                                   aahd_sweep_chunk.IsSurfaceSourceActive());
  async_comm_.WaitForUpstreamPsi();

  aahd_fluds->CopyNonLocalIncomingPsiToDevice();
  aahd_fluds->AllocateSaveAngularFlux(aahd_sweep_chunk.GetProblem(),
                                      aahd_sweep_chunk.GetGroupset());
  aahd_sweep_chunk.Sweep(*this);
  aahd_fluds->CopyPsiFromDevice();
  aahd_fluds->CopySaveAngularFluxFromDevice();
  stream_.synchronize();

  async_comm_.SendDownstreamPsi(static_cast<int>(this->GetID()));
  if (has_reflecting_boundaries_)
    aahd_fluds->CopyBoundaryPsiToAngleSet(aahd_sweep_chunk.GetGrid(), *this);
  for (auto& following_as : following_angle_sets_)
    following_as->starting_latch_->count_down();

  aahd_fluds->CopySaveAngularFluxToDestinationPsi(
    aahd_sweep_chunk.GetProblem(), aahd_sweep_chunk.GetGroupset(), *this);
  async_comm_.WaitForDownstreamPsi();

  async_comm_.WaitForDelayedIncomingPsi();

  executed_ = true;
  return AngleSetStatus::FINISHED;
}

} // namespace opensn
