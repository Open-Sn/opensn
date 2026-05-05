// SPDX-FileCopyrightText: 2026 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep/angle_set/cbcd_angle_set.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep/communicators/cbcd_async_comm.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep/fluds/cbcd_fluds.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep/spds/cbc.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep_chunks/sweep_chunk.h"
#include "framework/mesh/mesh_continuum/mesh_continuum.h"
#include "framework/data_types/range.h"
#include "framework/logging/log.h"
#include "framework/runtime.h"
#include "caliper/cali.h"

namespace opensn
{

CBCD_AngleSet::CBCD_AngleSet(size_t id,
                             size_t num_groups,
                             const SPDS& spds,
                             std::shared_ptr<FLUDS>& fluds,
                             const std::vector<size_t>& angle_indices,
                             std::map<uint64_t, std::shared_ptr<SweepBoundary>>& boundaries,
                             const MPICommunicatorSet& comm_set)
  : AngleSet(id, num_groups, spds, fluds, angle_indices, boundaries),
    cbc_spds_(dynamic_cast<const CBC_SPDS&>(spds)),
    async_comm_(id, *fluds, comm_set),
    stream_(),
    device_angle_indices_(angles_.size())
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
}

AsynchronousCommunicator*
CBCD_AngleSet::GetCommunicator()
{
  return static_cast<AsynchronousCommunicator*>(&async_comm_);
}

AngleSetStatus
CBCD_AngleSet::AngleSetAdvance(SweepChunk& sweep_chunk, AngleSetStatus permission)
{
  OpenSnLogicalError("CBCD_AngleSet::AngleSetAdvance should not be called. Routine is handled by "
                     "SweepScheduler::ScheduleAlgoAsyncFIFO.");
}

void
CBCD_AngleSet::ResetSweepBuffers()
{
  current_task_list_.clear();
  async_comm_.Reset();
  fluds_->ClearLocalAndReceivePsi();
  executed_ = false;
}

const double*
CBCD_AngleSet::PsiBoundary(uint64_t boundary_id,
                           unsigned int angle_num,
                           uint64_t cell_local_id,
                           unsigned int face_num,
                           unsigned int fi,
                           unsigned int g,
                           bool surface_source_active)
{
  if (boundaries_[boundary_id]->IsReflecting())
    return boundaries_[boundary_id]->PsiIncoming(cell_local_id, face_num, fi, angle_num, g);

  if (not surface_source_active)
    return boundaries_[boundary_id]->ZeroFlux(g);

  return boundaries_[boundary_id]->PsiIncoming(cell_local_id, face_num, fi, angle_num, g);
}

double*
CBCD_AngleSet::PsiReflected(uint64_t boundary_id,
                            unsigned int angle_num,
                            uint64_t cell_local_id,
                            unsigned int face_num,
                            unsigned int fi)
{
  return boundaries_[boundary_id]->PsiOutgoing(cell_local_id, face_num, fi, angle_num);
}

} // namespace opensn
