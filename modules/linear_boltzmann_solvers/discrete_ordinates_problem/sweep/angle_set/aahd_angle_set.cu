// SPDX-FileCopyrightText: 2026 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep/angle_set/aahd_angle_set.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep_chunks/aahd_sweep_chunk.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep/fluds/aahd_fluds.h"
#include "caliper/cali.h"
#include <algorithm>

namespace opensn
{

std::mutex AAHD_AngleSet::m;

AAHD_AngleSet::AAHD_AngleSet(size_t id,
                             const LBSGroupset& groupset,
                             const SPDS& spds,
                             std::shared_ptr<FLUDS>& fluds,
                             std::vector<size_t>& angle_indices,
                             std::map<uint64_t, std::shared_ptr<SweepBoundary>>& boundaries,
                             int maximum_message_size,
                             const MPICommunicatorSet& comm_set)
  : AngleSet(id, groupset, spds, fluds, angle_indices, boundaries),
    async_comm_(*fluds, num_groups_, angle_indices.size(), maximum_message_size, comm_set),
    device_angle_indices_(angles_.size())
{
  stream_ = crb::Stream();
  auto aahd_fluds = std::dynamic_pointer_cast<AAHD_FLUDS>(fluds_);
  aahd_fluds->GetStream() = stream_;
  aahd_fluds->GetCommonData().AddAssociatedAngleSet(this);
  crb::MemoryPinningManager<std::uint32_t> pin(angles_);
  crb::copy(device_angle_indices_, pin, angles_.size());
}

void
AAHD_AngleSet::InitializeDelayedUpstreamData()
{
  auto* aahd_fluds = static_cast<AAHD_FLUDS*>(fluds_.get());
  aahd_fluds->AllocateDelayedPrelocIOutgoingPsi();
  aahd_fluds->AllocateDelayedLocalPsi();
}

void
AAHD_AngleSet::PrepostReceives()
{
  async_comm_.PrepostReceiveUpstreamPsi(static_cast<int>(this->GetID()));
  async_comm_.PrepostReceiveDelayedData(static_cast<int>(this->GetID()));
  auto* aahd_fluds = static_cast<AAHD_FLUDS*>(fluds_.get());
  aahd_fluds->CopyDelayedPsiToDevice();
}

void
AAHD_AngleSet::WaitForDownstreamAndDelayed()
{
  async_comm_.WaitForDownstreamPsi();
  async_comm_.WaitForDelayedIncomingPsi();
}

bool
AAHD_AngleSet::IsReady()
{
  return (async_comm_.TestReceiveUpstreamPsi()) && (dependency_counter_ == 0);
}

AngleSetStatus
AAHD_AngleSet::AngleSetAdvance(SweepChunk& sweep_chunk, AngleSetStatus permission)
{
  CALI_CXX_MARK_SCOPE("AAHD_AngleSet::AngleSetAdvance");

  if (executed_)
    return AngleSetStatus::FINISHED;

  auto* aahd_fluds = static_cast<AAHD_FLUDS*>(fluds_.get());
  auto& aahd_sweep_chunk = static_cast<AAHDSweepChunk&>(sweep_chunk);

  aahd_fluds->CopyNonLocalIncomingPsiToDevice();
  aahd_fluds->AllocateSaveAngularFlux(aahd_sweep_chunk.GetProblem(),
                                      aahd_sweep_chunk.GetGroupset());
  aahd_sweep_chunk.Sweep(*this);
  aahd_fluds->CopyPsiFromDevice();
  aahd_fluds->CopySaveAngularFluxFromDevice();
  stream_.synchronize();

  if (not following_angle_sets_.empty())
  {
    std::scoped_lock lk(m);
    for (auto& following_as : following_angle_sets_)
      following_as->DecrementCounter();
  }
  async_comm_.SendDownstreamPsi(static_cast<int>(this->GetID()));

  aahd_fluds->CopySaveAngularFluxToDestinationPsi(
    aahd_sweep_chunk.GetProblem(), aahd_sweep_chunk.GetGroupset(), *this);

  executed_ = true;
  return AngleSetStatus::FINISHED;
}

void
AAHD_AngleSet::SyncDeviceAngleIndices()
{
  crb::MemoryPinningManager<std::uint32_t> pin(angles_);
  crb::copy(device_angle_indices_, pin, angles_.size());
}

void
AAHD_AngleSet::CopyBoundaryOffsetToDevice()
{
  host_boundary_offsets_.shrink_to_fit();
  if (host_boundary_offsets_.empty())
    return;
  crb::MemoryPinningManager<std::uint64_t> pin(host_boundary_offsets_);
  device_boundary_offsets_ = crb::DeviceMemory<std::uint64_t>(host_boundary_offsets_.size());
  crb::copy(device_boundary_offsets_, pin, host_boundary_offsets_.size());
}

} // namespace opensn
