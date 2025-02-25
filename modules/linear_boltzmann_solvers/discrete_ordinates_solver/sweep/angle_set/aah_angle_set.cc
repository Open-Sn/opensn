// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "modules/linear_boltzmann_solvers/discrete_ordinates_solver/sweep/angle_set/aah_angle_set.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_solver/sweep_chunks/sweep_chunk.h"
#include "framework/logging/log.h"
#include "framework/runtime.h"
#include "caliper/cali.h"

namespace opensn
{

AAH_AngleSet::AAH_AngleSet(size_t id,
                           size_t num_groups,
                           const SPDS& spds,
                           std::shared_ptr<FLUDS>& fluds,
                           std::vector<size_t>& angle_indices,
                           std::map<uint64_t, std::shared_ptr<SweepBoundary>>& boundaries,
                           int maximum_message_size,
                           const MPICommunicatorSet& comm_set)
  : AngleSet(id, num_groups, spds, fluds, angle_indices, boundaries),
    async_comm_(*fluds, num_groups_, angle_indices.size(), maximum_message_size, comm_set)
{
}

void
AAH_AngleSet::InitializeDelayedUpstreamData()
{
  async_comm_.InitializeDelayedUpstreamData();
}

AngleSetStatus
AAH_AngleSet::AngleSetAdvance(SweepChunk& sweep_chunk, AngleSetStatus permission)
{
  CALI_CXX_MARK_SCOPE("AAH_AngleSet::AngleSetAdvance");

  if (executed_)
  {
    if (not async_comm_.DoneSending())
      async_comm_.ClearDownstreamBuffers();
    return AngleSetStatus::FINISHED;
  }

  // Check upstream data available
  AngleSetStatus status = async_comm_.ReceiveUpstreamPsi(static_cast<int>(this->GetID()));

  // Also check boundaries
  for (auto& [bid, boundary] : boundaries_)
    if (not boundary->CheckAnglesReadyStatus(angles_))
    {
      status = AngleSetStatus::RECEIVING;
      break;
    }

  if (status == AngleSetStatus::RECEIVING)
    return status;
  else if (status == AngleSetStatus::READY_TO_EXECUTE and permission == AngleSetStatus::EXECUTE)
  {
    async_comm_.InitializeLocalAndDownstreamBuffers();

    sweep_chunk.Sweep(*this); // Execute chunk

    // Send outgoing psi and clear local and receive buffers
    async_comm_.SendDownstreamPsi(static_cast<int>(this->GetID()));
    async_comm_.ClearLocalAndReceiveBuffers();

    // Update boundary readiness
    for (auto& [bid, boundary] : boundaries_)
      boundary->UpdateAnglesReadyStatus(angles_);

    executed_ = true;
    return AngleSetStatus::FINISHED;
  }
  else
    return AngleSetStatus::READY_TO_EXECUTE;
}

AngleSetStatus
AAH_AngleSet::FlushSendBuffers()
{
  if (not async_comm_.DoneSending())
    async_comm_.ClearDownstreamBuffers();

  if (async_comm_.DoneSending())
    return AngleSetStatus::MESSAGES_SENT;

  return AngleSetStatus::MESSAGES_PENDING;
}

int
AAH_AngleSet::GetMaxBufferMessages() const
{
  return async_comm_.GetMaxNumMessages();
}

void
AAH_AngleSet::SetMaxBufferMessages(int count)
{
  async_comm_.SetMaxNumMessages(count);
}

void
AAH_AngleSet::ResetSweepBuffers()
{
  async_comm_.Reset();
  executed_ = false;
}

bool
AAH_AngleSet::ReceiveDelayedData()
{
  return async_comm_.ReceiveDelayedData(static_cast<int>(this->GetID()));
}

const double*
AAH_AngleSet::PsiBoundary(uint64_t boundary_id,
                          unsigned int angle_num,
                          uint64_t cell_local_id,
                          unsigned int face_num,
                          unsigned int fi,
                          int g,
                          bool surface_source_active)
{
  if (boundaries_[boundary_id]->IsReflecting())
    return boundaries_[boundary_id]->PsiIncoming(cell_local_id, face_num, fi, angle_num, g);

  if (not surface_source_active)
    return boundaries_[boundary_id]->ZeroFlux(g);

  return boundaries_[boundary_id]->PsiIncoming(cell_local_id, face_num, fi, angle_num, g);
}

double*
AAH_AngleSet::PsiReflected(uint64_t boundary_id,
                           unsigned int angle_num,
                           uint64_t cell_local_id,
                           unsigned int face_num,
                           unsigned int fi)
{
  return boundaries_[boundary_id]->PsiOutgoing(cell_local_id, face_num, fi, angle_num);
}

} // namespace opensn
