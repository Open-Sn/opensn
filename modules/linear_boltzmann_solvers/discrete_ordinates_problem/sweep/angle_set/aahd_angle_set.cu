// SPDX-FileCopyrightText: 2026 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep/angle_set/aahd_angle_set.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep_chunks/aahd_sweep_chunk.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep/fluds/aahd_fluds.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep/spds/aah.h"
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
                             const SweepCommunicator& sweep_communicator)
  : AngleSet(id, groupset, spds, fluds, angle_indices, boundaries),
    async_comm_(*fluds,
                groupset.id,
                id,
                num_groups_,
                angle_indices.size(),
                maximum_message_size,
                sweep_communicator),
    device_angle_indices_(angles_.size())
{
  stream_ = crb::Stream();
  auto aahd_fluds = std::dynamic_pointer_cast<AAHD_FLUDS>(fluds_);
  aahd_fluds->GetStream() = stream_;
  aahd_fluds->GetCommonData().AddAssociatedAngleSet(this);
  AAHD_AngleSet::SyncDeviceAngleIndices();
}

void
AAHD_AngleSet::InitializeDelayedUpstreamData()
{
  auto* aahd_fluds = dynamic_cast<AAHD_FLUDS*>(fluds_.get());
  aahd_fluds->AllocateDelayedPrelocIOutgoingPsi();
  aahd_fluds->AllocateDelayedLocalPsi();
}

void
AAHD_AngleSet::PrepostReceives()
{
  async_comm_.PrepostReceiveUpstreamPsi();
  async_comm_.PrepostReceiveDelayedData();
  auto* aahd_fluds = dynamic_cast<AAHD_FLUDS*>(fluds_.get());
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
  if (executed_)
    return AngleSetStatus::FINISHED;

  auto* aahd_fluds = dynamic_cast<AAHD_FLUDS*>(fluds_.get());
  auto& aahd_sweep_chunk = dynamic_cast<AAHDSweepChunk&>(sweep_chunk);

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
  async_comm_.SendDownstreamPsi();

  aahd_fluds->CopySaveAngularFluxToDestinationPsi(
    aahd_sweep_chunk.GetProblem(), aahd_sweep_chunk.GetGroupset(), *this);

  executed_ = true;
  return AngleSetStatus::FINISHED;
}

void
AAHD_AngleSet::SyncDeviceAngleIndices()
{
  crb::copy(device_angle_indices_, angles_, angles_.size());
}

void
AAHD_AngleSet::CopyBoundaryOffsetToDevice()
{
  if (host_boundary_offsets_.empty())
  {
    device_boundary_offsets_.reset();
    return;
  }
  device_boundary_offsets_ = crb::DeviceMemory<std::uint64_t>(host_boundary_offsets_.size());
  crb::copy(device_boundary_offsets_, host_boundary_offsets_, host_boundary_offsets_.size());
  host_boundary_offsets_.clear();
}

void
AAHD_AngleSet::SetFLUDSStream(crb::Stream& stream)
{
  dynamic_cast<AAHD_FLUDS*>(fluds_.get())->GetStream() = stream;
}

void
AAHD_AngleSet::PrepareForBatch(SweepChunk& sweep_chunk, crb::Stream& stream)
{
  auto* aahd_fluds = dynamic_cast<AAHD_FLUDS*>(fluds_.get());
  auto& aahd_sweep_chunk = dynamic_cast<AAHDSweepChunk&>(sweep_chunk);
  // Put this angle set's FLUDS copies on its run's lane stream, so they stay ordered relative to
  // that run's fused kernel launches.
  aahd_fluds->GetStream() = stream;
  aahd_fluds->CopyNonLocalIncomingPsiToDevice();
  aahd_fluds->AllocateSaveAngularFlux(aahd_sweep_chunk.GetProblem(),
                                      aahd_sweep_chunk.GetGroupset());
}

std::size_t
AAHD_AngleSet::NumLevels() const
{
  const auto& spds = dynamic_cast<const AAH_SPDS&>(spds_);
  return spds.GetLevelizedLocalSubgrid().size();
}

void
AAHD_AngleSet::FillBatchSlots(std::size_t num_levels,
                              std::size_t slot_stride,
                              gpu_kernel::AAHD_BatchSlot* slots) const
{
  auto* aahd_fluds = dynamic_cast<AAHD_FLUDS*>(fluds_.get());
  const auto& spds = dynamic_cast<const AAH_SPDS&>(spds_);
  const auto& levelized_spls = spds.GetLevelizedLocalSubgrid();

  gpu_kernel::AAHD_BatchSlot slot;
  slot.directions = device_angle_indices_.get();
  slot.boundary_offset = device_boundary_offsets_.get();
  slot.angleset_size = static_cast<std::uint32_t>(angles_.size());
  slot.flud_index = aahd_fluds->GetCommonData().GetDeviceIndex();
  slot.flud_data = aahd_fluds->GetDevicePointerSet();
  slot.saved_psi = aahd_fluds->GetSavedAngularFluxDevicePointer();
  for (std::size_t level = 0; level < num_levels; ++level)
  {
    if (level < levelized_spls.size())
    {
      slot.cells_to_sweep = spds.GetDeviceLevelVector(level);
      slot.num_cells = static_cast<unsigned int>(levelized_spls[level].size());
    }
    else
    {
      slot.cells_to_sweep = nullptr;
      slot.num_cells = 0;
    }
    slots[level * slot_stride] = slot;
  }
}

void
AAHD_AngleSet::SubmitCompletionCopies()
{
  auto* aahd_fluds = dynamic_cast<AAHD_FLUDS*>(fluds_.get());
  aahd_fluds->CopyPsiFromDevice();
  aahd_fluds->CopySaveAngularFluxFromDevice();
}

void
AAHD_AngleSet::FinishCompletion(SweepChunk& sweep_chunk)
{
  auto* aahd_fluds = dynamic_cast<AAHD_FLUDS*>(fluds_.get());
  auto& aahd_sweep_chunk = dynamic_cast<AAHDSweepChunk&>(sweep_chunk);

  // The caller saw this angle set's lane stream complete after SubmitCompletionCopies, so the
  // host buffers sent below and read by CopySaveAngularFluxToDestinationPsi are complete.

  if (not following_angle_sets_.empty())
  {
    std::scoped_lock lk(m);
    for (auto& following_as : following_angle_sets_)
      following_as->DecrementCounter();
  }
  async_comm_.SendDownstreamPsi();

  aahd_fluds->CopySaveAngularFluxToDestinationPsi(
    aahd_sweep_chunk.GetProblem(), aahd_sweep_chunk.GetGroupset(), *this);

  executed_ = true;
}

} // namespace opensn
