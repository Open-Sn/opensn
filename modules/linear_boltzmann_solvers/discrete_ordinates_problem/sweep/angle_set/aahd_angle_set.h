// SPDX-FileCopyrightText: 2026 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep/angle_set/angle_set.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep/communicators/aahd_async_comm.h"
#include "modules/linear_boltzmann_solvers/lbs_problem/device/memory_pinner.h"
#include "caribou/main.hpp"
#include <latch>
#include <memory>
#include <set>

namespace crb = caribou;

namespace opensn
{

/// AAHD angle set.
class AAHD_AngleSet : public AngleSet
{
public:
  AAHD_AngleSet(size_t id,
                size_t num_groups,
                const SPDS& spds,
                std::shared_ptr<FLUDS>& fluds,
                std::vector<size_t>& angle_indices,
                std::map<uint64_t, std::shared_ptr<SweepBoundary>>& boundaries,
                int maximum_message_size,
                const MPICommunicatorSet& in_comm_set);

  crb::Stream& GetStream() { return stream_; }

  void InitializeDelayedUpstreamData() override;

  int GetMaxBufferMessages() const override { return async_comm_.GetMaxNumMessages(); }

  void SetMaxBufferMessages(int count) override { async_comm_.SetMaxNumMessages(count); }

  AngleSetStatus AngleSetAdvance(SweepChunk& sweep_chunk, AngleSetStatus permission) override;

  AngleSetStatus FlushSendBuffers() override { return AngleSetStatus::MESSAGES_SENT; }

  void ResetSweepBuffers() override;

  bool ReceiveDelayedData() override { return true; }

  const double* PsiBoundary(uint64_t boundary_id,
                            unsigned int angle_num,
                            uint64_t cell_local_id,
                            unsigned int face_num,
                            unsigned int fi,
                            unsigned int g,
                            bool surface_source_active) override
  {
    if (boundaries_[boundary_id]->IsReflecting())
      return boundaries_[boundary_id]->PsiIncoming(cell_local_id, face_num, fi, angle_num, g);

    if (not surface_source_active)
      return boundaries_[boundary_id]->ZeroFlux(g);

    return boundaries_[boundary_id]->PsiIncoming(cell_local_id, face_num, fi, angle_num, g);
  }

  double* PsiReflected(uint64_t boundary_id,
                       unsigned int angle_num,
                       uint64_t cell_local_id,
                       unsigned int face_num,
                       unsigned int fi) override
  {
    return boundaries_[boundary_id]->PsiOutgoing(cell_local_id, face_num, fi, angle_num);
  }

protected:
  /// Associated CUDA stream.
  crb::Stream stream_;

  /// Asynchronous communicator.
  AAHD_ASynchronousCommunicator async_comm_;
};

} // namespace opensn
