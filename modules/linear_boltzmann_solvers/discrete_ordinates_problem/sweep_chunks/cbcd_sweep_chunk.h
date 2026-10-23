// SPDX-FileCopyrightText: 2026 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep/angle_set/cbcd_angle_set.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep/communicators/cbcd_async_comm.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep/fluds/cbcd_fluds.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep_chunks/sweep_chunk.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/discrete_ordinates_problem.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep_chunks/gpu_kernel/arguments.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep_chunks/gpu_kernel/cbcd_batch_slot.h"
#include "caribou/main.hpp"

namespace crb = caribou;

namespace opensn
{

/// Device CBC sweep chunk.
class CBCDSweepChunk : public SweepChunk
{
public:
  CBCDSweepChunk(DiscreteOrdinatesProblem& problem, LBSGroupset& groupset);

  /// Stop the communicator before the angle-set storage it reads is destroyed.
  ~CBCDSweepChunk() override;

  DiscreteOrdinatesProblem& GetProblem() const { return problem_; }

  const LBSGroupset& GetGroupset() const { return groupset_; }

  const std::vector<CBCD_AngleSet*>& GetAngleSets() const { return angle_sets_; }

  /// Return true when lagged boundary angular fluxes must be uploaded before each sweep.
  bool NeedsBoundaryUpload() const { return needs_boundary_upload_; }

  /**
   * Start a sweep on the communicator's progress thread.
   * This also creates the lanes on the first call or when the number of workers changes.
   *
   * \param num_workers Number of sweep workers.
   */
  void StartCommunicator(std::size_t num_workers);

  /// Wait until the communicator has sent and received all data of the sweep.
  void StopCommunicator();

  /// Rebuild the sweep arguments of every angle set and copy them to the device.
  void RefreshKernelArguments();

  /**
   * \name Fused launches
   * Each sweep worker has one lane: a stream with at most one fused launch in flight. During a
   * worker visit, the worker's angle sets add their ready cells as batch slots to the lane's next
   * launch, which SubmitFusedLaunch launches at the end of the visit.
   * \{
   */
  /**
   * Observe whether the worker's launch in flight has completed.
   * Call this once at the start of each worker visit. A completed launch's batch slots are
   * released, and IsLaneLaunchCompleted returns true until the next call.
   *
   * \param worker_id Sweep worker index.
   */
  void PollLane(std::size_t worker_id) { lanes_[worker_id]->Poll(); }

  /// Return true when the last PollLane of the worker observed the completion of its launch.
  bool IsLaneLaunchCompleted(std::size_t worker_id) const { return lanes_[worker_id]->completed; }

  /// Return true while the worker has a launch in flight.
  bool IsLaneInFlight(std::size_t worker_id) const { return lanes_[worker_id]->in_flight; }

  /**
   * Add an angle set's ready cells to its worker's next launch.
   *
   * \param worker_id Sweep worker that owns the angle set.
   * \param angle_set_id Angle-set index.
   * \param local_cell_ids Local IDs of the ready cells, which is readable by the device until the
   * launch completes.
   * \param num_ready_cells Number of ready cells.
   */
  void AddBatchSlot(std::size_t worker_id,
                    std::size_t angle_set_id,
                    const std::uint32_t* local_cell_ids,
                    std::uint32_t num_ready_cells);

  /**
   * Launch the batch slots added to a worker's lane since its previous launch.
   * Batches are not launched while the lane has a launch in flight or has no batch slots.
   *
   * \param worker_id Sweep worker index.
   */
  void SubmitFusedLaunch(std::size_t worker_id);
  /** \} */

private:
  /**
   * Stream and batch slots of one sweep worker's fused launches.
   * The lane is neither copyable nor destroyed before its stream completes
   * because the batch slots must not change while a launch is in flight.
   */
  struct FusedLane
  {
    FusedLane() = default;
    FusedLane(const FusedLane&) = delete;
    FusedLane& operator=(const FusedLane&) = delete;
    ~FusedLane() { stream.synchronize(); }

    /// Stream of the worker's launches.
    crb::Stream stream;
    /// Batch slots of the next or in-flight launch, which is read by the device.
    crb::MappedHostVector<gpu_kernel::CBCD_BatchSlot> batch_slots;
    /// True from submission until Poll observes that the launch has completed.
    bool in_flight = false;
    /// True during the worker visit in which Poll observed completion of the launch.
    bool completed = false;

    /// Observe completion of the in-flight launch and release its batch slots.
    void Poll()
    {
      completed = in_flight and stream.is_completed();
      if (completed)
      {
        in_flight = false;
        batch_slots.clear();
      }
    }
  };

  /// Launch geometry and buffers of one angle set.
  struct AngleSetLaunch
  {
    /// Threads per block along the angle-group (x) and cell (y) axes.
    crb::Dim3 threads_per_block;
    /// Number of blocks spanning the angle-group axis.
    unsigned int num_stride_blocks = 0;
    /// FLUDS of the angle set.
    CBCD_FLUDS* fluds = nullptr;
    /// Saved-angular-flux device buffer of the angle set (nullptr if fluxes are not saved).
    double* device_saved_psi = nullptr;
  };

  DiscreteOrdinatesProblem& problem_;
  /// True when a boundary has lagged angular fluxes for this groupset.
  bool needs_boundary_upload_ = false;
  /// Communicator shared by every angle set of the groupset.
  std::unique_ptr<CBCD_AsynchronousCommunicator> async_comm_;
  /// Angle sets of the groupset.
  std::vector<CBCD_AngleSet*> angle_sets_;
  /// Launch data indexed by angle set.
  std::vector<AngleSetLaunch> angle_set_launches_;
  /// Host copy of the sweep arguments indexed by angle set.
  crb::HostVector<gpu_kernel::Arguments<SweepKind::CBC>> host_arguments_;
  /// Device copy of the sweep arguments indexed by angle set.
  crb::DeviceMemory<gpu_kernel::Arguments<SweepKind::CBC>> device_arguments_;
  /// Lanes indexed by sweep worker.
  std::vector<std::unique_ptr<FusedLane>> lanes_;
};

} // namespace opensn
