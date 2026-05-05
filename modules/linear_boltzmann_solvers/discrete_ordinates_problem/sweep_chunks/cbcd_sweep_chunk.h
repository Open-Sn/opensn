// SPDX-FileCopyrightText: 2026 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep/angle_set/cbcd_angle_set.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep/fluds/cbcd_fluds.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep_chunks/sweep_chunk.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/discrete_ordinates_problem.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep_chunks/gpu_kernel/arguments.h"
#include "caribou/main.hpp"

namespace crb = caribou;

namespace opensn
{

/// CBC sweep chunk for device.
class CBCDSweepChunk : public SweepChunk
{
public:
  CBCDSweepChunk(DiscreteOrdinatesProblem& problem, LBSGroupset& groupset);

  DiscreteOrdinatesProblem& GetProblem() const { return problem_; }

  const LBSGroupset& GetGroupset() const { return groupset_; }

  unsigned int GetGroupsetGroupIndex() const { return groupset_.first_group; }

  const CellLBSView& GetCellTransportView(std::uint64_t cell_local_id) const
  {
    return cell_transport_views_[cell_local_id];
  }

  const std::vector<CBCD_AngleSet*>& GetAngleSets() const { return angle_sets_; }

  const std::vector<CBCD_FLUDS*>& GetFLUDS() const { return fluds_list_; }

  std::vector<crb::Stream>& GetStreams() { return streams_list_; }

  using SweepChunk::Sweep;
  void Sweep(const std::vector<std::uint32_t>& cell_local_ids, size_t angle_set_id);

private:
  DiscreteOrdinatesProblem& problem_;
  std::vector<CBCD_AngleSet*> angle_sets_;
  std::vector<CBCD_FLUDS*> fluds_list_;
  std::vector<crb::Stream> streams_list_;
  std::vector<gpu_kernel::Arguments<gpu_kernel::SweepType::CBC>> kernel_args_list_;
  std::vector<crb::Dim3> block_sizes_;
  std::vector<unsigned int> grid_size_x_list_;
};

} // namespace opensn
