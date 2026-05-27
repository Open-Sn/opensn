// SPDX-FileCopyrightText: 2026 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep/angle_set/cbcd_angle_set.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep/communicators/cbcd_async_comm.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep/fluds/cbcd_fluds.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep_chunks/sweep_chunk.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/discrete_ordinates_problem.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep_chunks/gpu_kernel/arguments.h"
#include "caribou/main.hpp"

namespace crb = caribou;

namespace opensn
{

/// CBCD sweep chunk.
class CBCDSweepChunk : public SweepChunk
{
public:
  /// Build persistent kernel launches and the groupset-wide communicator.
  CBCDSweepChunk(DiscreteOrdinatesProblem& problem, LBSGroupset& groupset);

  /// Stop the communicator before destroying angle-set storage.
  ~CBCDSweepChunk() override;

  /// Return the owning transport problem.
  DiscreteOrdinatesProblem& GetProblem() const { return problem_; }

  /// Return the active groupset.
  const LBSGroupset& GetGroupset() const { return groupset_; }

  /// Return the first global group index in the active groupset.
  unsigned int GetGroupsetGroupIndex() const { return groupset_.first_group; }

  /// Return transport metadata for one local cell.
  const CellLBSView& GetCellTransportView(std::uint64_t cell_local_id) const
  {
    return cell_transport_views_[cell_local_id];
  }

  /// Return CBCD angle sets in scheduler order.
  const std::vector<CBCD_AngleSet*>& GetAngleSets() const { return angle_sets_; }

  /// Start the MPI progress thread and configure worker-owned queues.
  void StartCommunicator(std::size_t num_workers);

  /// Drain and stop the MPI progress thread.
  void StopCommunicator();

  /// Refresh problem-dependent arguments cached for each angle set.
  void RefreshKernelArguments();

  using SweepChunk::Sweep;
  /// Launch one ready-cell batch.
  void Sweep(std::uint32_t num_ready_cells,
             std::size_t angle_set_id,
             const std::uint32_t* local_cell_ids);

private:
  struct KernelLaunch
  {
    /// Persistent kernel arguments refreshed when problem vectors change.
    gpu_kernel::Arguments<SweepKind::CBC> arguments;
    /// Fixed block geometry and stride-axis grid extent.
    crb::Dim3 threads_per_block;
    unsigned int num_stride_blocks;
    /// Owning FLUDS and optional saved-psi device storage.
    CBCD_FLUDS* fluds;
    double* device_saved_psi;
  };
  /// Owning problem and groupset-wide aggregated communicator.
  DiscreteOrdinatesProblem& problem_;
  std::unique_ptr<CBCD_AsynchronousCommunicator> async_comm_;
  /// Angle sets and their persistent kernel launches in scheduler order.
  std::vector<CBCD_AngleSet*> angle_sets_;
  std::vector<KernelLaunch> kernel_launches_;
};

} // namespace opensn
