// SPDX-FileCopyrightText: 2026 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep/angle_set/aahd_angle_set.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep_chunks/aahd_sweep_chunk.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep/fluds/aahd_fluds.h"
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
  std::dynamic_pointer_cast<AAHD_FLUDS>(fluds_)->GetStream() = stream_;
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

  // Note: Order matters. Simultaneous operations (thread, stream, mpi) must not be separated by
  // empty lines. Empty lines represent synchronization points with thread.

  // stream: allocate memory and asynchronously copy delayed psi to device
  aahd_fluds->AllocateInternalLocalPsi();
  aahd_fluds->AllocateOutgoingPsi();
  aahd_fluds->CopyDelayedPsiToDevice();
  // mpi: allocate memory and pre-post receive for non-local incoming psi
  aahd_fluds->AllocatePrelocIOutgoingPsi();
  async_comm_.PrepostReceiveUpstreamPsi(static_cast<int>(this->GetID()));
  // thread: wait for reflecting boundary data and copy boundary data to device
  starting_latch_->wait();
  aahd_fluds->CopyBoundaryToDevice(aahd_sweep_chunk.GetGrid(),
                                   *this,
                                   aahd_sweep_chunk.GetGroupset(),
                                   aahd_sweep_chunk.IsSurfaceSourceActive());

  // thread: wait for upstream data to arrive
  async_comm_.WaitForUpstreamPsi();

  // stream: copy data to device, sweep and copy data back to host buffers
  aahd_fluds->CopyNonLocalIncomingPsiToDevice();
  aahd_fluds->AllocateSaveAngularFlux(aahd_sweep_chunk.GetProblem(),
                                      aahd_sweep_chunk.GetGroupset());
  sweep_chunk.Sweep(*this);
  aahd_fluds->CopyPsiFromDevice();
  // thread: wait for all sweep computations and data transfer to finish
  stream_.synchronize();

  // stream: copy save angular flux to contiguous buffer on host, free local + upstream
  aahd_fluds->CopySaveAngularFluxFromDevice();
  aahd_fluds->ClearLocalAndReceivePsi();
  // mpi: send non-local outgoing psi and pre-post receive delayed data
  async_comm_.SendDownstreamPsi(static_cast<int>(this->GetID()));
  async_comm_.PrepostReceiveDelayedData(static_cast<int>(this->GetID()));
  // thread: copy boundary data to angle set to unlock the latches of following angle sets,
  aahd_fluds->CopyBoundaryPsiToAngleSet(aahd_sweep_chunk.GetGrid(), *this);
  for (auto& following_as : following_angle_sets_)
  {
    following_as->starting_latch_->count_down();
  }
  // thread: copy save angular flux to destination psi,
  if (aahd_fluds->HasSaveAngularFlux())
  {
    stream_.synchronize();
    aahd_fluds->CopySaveAngularFluxToDestinationPsi(
      aahd_sweep_chunk.GetProblem(), aahd_sweep_chunk.GetGroupset(), *this);
  }
  // thread: wait for downstream sends
  async_comm_.WaitForDownstreamPsi();

  // stream: deallocate downstream psi memory
  aahd_fluds->ClearSendPsi();
  // thread: wait for delayed incoming messages
  async_comm_.WaitForDelayedIncomingPsi();
  // thread: wait for all device memory deallocations
  stream_.synchronize();

  executed_ = true;
  return AngleSetStatus::FINISHED;
}

} // namespace opensn
