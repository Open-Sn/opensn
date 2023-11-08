#include "framework/mesh/SweepUtilities/AngleSet/AAH_AngleSet.h"

#include "framework/mesh/SweepUtilities/sweepchunk_base.h"

#include "framework/chi_runtime.h"
#include "framework/logging/chi_log.h"

namespace chi_mesh::sweep_management
{

AAH_AngleSet::AAH_AngleSet(size_t id,
                           size_t in_numgrps,
                           size_t in_ref_subset,
                           const SPDS& in_spds,
                           std::shared_ptr<FLUDS>& in_fluds,
                           std::vector<size_t>& angle_indices,
                           std::map<uint64_t, std::shared_ptr<SweepBndry>>& sim_boundaries,
                           int sweep_eager_limit,
                           const chi::ChiMPICommunicatorSet& in_comm_set)
  : AngleSet(id, in_numgrps, in_spds, in_fluds, angle_indices, sim_boundaries, in_ref_subset),
    async_comm_(*in_fluds, num_grps, angle_indices.size(), sweep_eager_limit, in_comm_set)
{
}

void
AAH_AngleSet::InitializeDelayedUpstreamData()
{
  async_comm_.InitializeDelayedUpstreamData();
}

AngleSetStatus
AAH_AngleSet::AngleSetAdvance(SweepChunk& sweep_chunk,
                              const std::vector<size_t>& timing_tags,
                              ExecutionPermission permission)
{
  typedef AngleSetStatus Status;

  if (executed_)
  {
    if (!async_comm_.DoneSending()) async_comm_.ClearDownstreamBuffers();
    return AngleSetStatus::FINISHED;
  }

  // Check upstream data available
  Status status = async_comm_.ReceiveUpstreamPsi(static_cast<int>(this->GetID()));

  // Also check boundaries
  for (auto& [bid, bndry] : ref_boundaries_)
    if (not bndry->CheckAnglesReadyStatus(angles_, ref_group_subset_))
    {
      status = Status::RECEIVING;
      break;
    }

  if (status == Status::RECEIVING) return status;
  else if (status == Status::READY_TO_EXECUTE and permission == ExecutionPermission::EXECUTE)
  {
    async_comm_.InitializeLocalAndDownstreamBuffers();

    Chi::log.LogEvent(timing_tags[0], chi::ChiLog::EventType::EVENT_BEGIN);
    sweep_chunk.Sweep(*this); // Execute chunk
    Chi::log.LogEvent(timing_tags[0], chi::ChiLog::EventType::EVENT_END);

    // Send outgoing psi and clear local and receive buffers
    async_comm_.SendDownstreamPsi(static_cast<int>(this->GetID()));
    async_comm_.ClearLocalAndReceiveBuffers();

    // Update boundary readiness
    for (auto& [bid, bndry] : ref_boundaries_)
      bndry->UpdateAnglesReadyStatus(angles_, ref_group_subset_);

    executed_ = true;
    return AngleSetStatus::FINISHED;
  }
  else
    return AngleSetStatus::READY_TO_EXECUTE;
}

AngleSetStatus
AAH_AngleSet::FlushSendBuffers()
{
  if (!async_comm_.DoneSending()) async_comm_.ClearDownstreamBuffers();

  if (async_comm_.DoneSending()) return AngleSetStatus::MESSAGES_SENT;

  return AngleSetStatus::MESSAGES_PENDING;
}

int
AAH_AngleSet::GetMaxBufferMessages() const
{
  return async_comm_.max_num_mess;
}

void
AAH_AngleSet::SetMaxBufferMessages(int new_max)
{
  async_comm_.max_num_mess = new_max;
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
AAH_AngleSet::PsiBndry(uint64_t bndry_map,
                       unsigned int angle_num,
                       uint64_t cell_local_id,
                       unsigned int face_num,
                       unsigned int fi,
                       int g,
                       size_t gs_ss_begin,
                       bool surface_source_active)
{
  if (ref_boundaries_[bndry_map]->IsReflecting())
    return ref_boundaries_[bndry_map]->HeterogeneousPsiIncoming(
      cell_local_id, face_num, fi, angle_num, g, gs_ss_begin);

  if (not surface_source_active) return ref_boundaries_[bndry_map]->ZeroFlux(g);

  return ref_boundaries_[bndry_map]->HeterogeneousPsiIncoming(
    cell_local_id, face_num, fi, angle_num, g, gs_ss_begin);
}

double*
AAH_AngleSet::ReflectingPsiOutBoundBndry(uint64_t bndry_map,
                                         unsigned int angle_num,
                                         uint64_t cell_local_id,
                                         unsigned int face_num,
                                         unsigned int fi,
                                         size_t gs_ss_begin)
{
  return ref_boundaries_[bndry_map]->HeterogeneousPsiOutgoing(
    cell_local_id, face_num, fi, angle_num, gs_ss_begin);
}

} // namespace chi_mesh::sweep_management
