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

/**
 * CBCD sweep chunk.
 *
 * Owns the shared CBCD communicator for one groupset, caches per-angle-set kernel
 * launch parameters, and coordinates the transfer boundaries between the device sweep
 * kernels and the host-side CBCD scheduler.
 */
class CBCDSweepChunk : public SweepChunk
{
public:
  /**
   * Construct the CBCD sweep chunk for one groupset.
   *
   * \param problem Discrete ordinates problem owning the sweep state.
   * \param groupset Groupset served by this sweep chunk.
   */
  CBCDSweepChunk(DiscreteOrdinatesProblem& problem, LBSGroupset& groupset);

  ~CBCDSweepChunk() override;

  /// Return the discrete ordinates problem owning this sweep chunk.
  DiscreteOrdinatesProblem& GetProblem() const { return problem_; }

  /// Return the groupset served by this sweep chunk.
  const LBSGroupset& GetGroupset() const { return groupset_; }

  /// Return the first group index of the groupset.
  unsigned int GetGroupsetGroupIndex() const { return groupset_.first_group; }

  /// Return the cell transport view for one local cell.
  const CellLBSView& GetCellTransportView(std::uint64_t cell_local_id) const
  {
    return cell_transport_views_[cell_local_id];
  }

  /// Return the CBCD angle sets coordinated by this sweep chunk.
  const std::vector<CBCD_AngleSet*>& GetAngleSets() const { return angle_sets_; }

  /// Start the aggregated communicator thread.
  void StartCommunicator();

  /// Stop the aggregated communicator thread.
  void StopCommunicator();

  /// Refresh cached kernel arguments once at the start of a sweep.
  void RefreshCachedKernelArgs();

  using SweepChunk::Sweep;
  /**
   * Launch the CBC sweep kernel for one angle set.
   *
   * \param num_ready_cells Number of local cells in the batch.
   * \param angle_set_id Producing angle-set ID.
   * \param local_cell_ids Pointer to the mapped host cell-ID buffer for the batch.
   */
  void Sweep(std::uint32_t num_ready_cells,
             std::size_t angle_set_id,
             const std::uint32_t* local_cell_ids);

private:
  /// Cached launch data for one angle set.
  struct CachedKernelParams
  {
    /// Packed kernel arguments.
    gpu_kernel::Arguments<gpu_kernel::SweepType::CBC> args;
    /// Device block size for the launch.
    crb::Dim3 block_size;
    /// Device grid size in x.
    unsigned int grid_size_x;
    /// FLUDS instance bound to the angle set.
    CBCD_FLUDS* fluds;
    /// Device pointer to saved angular fluxes.
    double* device_saved_psi;
  };
  /// Owning reference to the discrete ordinates problem.
  DiscreteOrdinatesProblem& problem_;
  /// Aggregated communicator owned by this sweep chunk.
  std::unique_ptr<CBCD_AsynchronousCommunicator> async_comm_;
  /// Anglesets managed by this sweep chunk.
  std::vector<CBCD_AngleSet*> angle_sets_;
  /// Per-angleset cached kernel launch params.
  std::vector<CachedKernelParams> cached_params_;
};

} // namespace opensn
