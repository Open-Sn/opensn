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
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep_chunks/gpu_kernel/cbcd_sweep_launch.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep_chunks/gpu_kernel/round_up.h"
#include "caribou/main.hpp"
#include <array>

namespace crb = caribou;

namespace opensn
{

/// Device CBC sweep chunk for a groupset. Each worker owns one lane and its stream.
/// The problem, groupset, angle aggregation and FLUDS must outlive the chunk and its GPU/MPI work.
class CBCDSweepChunk : public SweepChunk
{
public:
  /// Construct from a built problem and groupset before sweep workers start.
  CBCDSweepChunk(DiscreteOrdinatesProblem& problem, LBSGroupset& groupset);

  DiscreteOrdinatesProblem& GetProblem() const { return problem_; }

  const LBSGroupset& GetGroupset() const { return groupset_; }

  const std::vector<CBCD_AngleSet*>& GetAngleSets() const { return angle_sets_; }

  /// Return true when lagged boundary angular fluxes must be uploaded before each sweep.
  bool NeedsBoundaryUpload() const { return needs_boundary_upload_; }

  std::size_t NumFusedLanes() const { return lanes_.size(); }

  void StartCommunicator();

  void StopCommunicator();

  /// Rebuild device sweep arguments on the caller thread, with all sweep workers idle.
  void RefreshKernelArguments();

  /**
   * \name Fused launches
   * Only the lane's worker may call its methods. A lane accepts new slots only while idle,
   * submits them on one stream at the end of a visit and polls the entire submission at the
   * start of each visit. worker_id must be less than NumFusedLanes(). Keep device arguments
   * and saved flux allocated through completion.
   * \{
   */
  /// Poll the lane and publish outgoing faces if the submission completes.
  void PollLane(std::size_t worker_id);

  /// Return true during the worker visit in which PollLane observed completion.
  bool IsLaneLaunchCompleted(std::size_t worker_id) const { return lanes_[worker_id]->completed; }

  bool IsLaneInFlight(std::size_t worker_id) const { return lanes_[worker_id]->in_flight; }

  /// Copy an angle set's ready cell IDs and their incoming payload offsets into the idle lane's
  /// next submission.
  void AddBatchSlot(std::size_t worker_id,
                    std::size_t angle_set_id,
                    const std::uint32_t* local_cell_ids,
                    std::uint32_t num_ready_cells);

  /// Submit the queued slots. Slots may span several kernels on the same stream.
  void SubmitFusedLaunch(std::size_t worker_id);
  /** \} */

private:
  /// Worker-owned stream and launch storage. Destroy only after the stream completes.
  struct FusedLane
  {
    FusedLane() = default;
    FusedLane(const FusedLane&) = delete;
    FusedLane& operator=(const FusedLane&) = delete;
    ~FusedLane() { stream.synchronize(); }

    crb::Stream stream;
#if defined(SYCL_LANGUAGE_VERSION) && defined(__INTEL_LLVM_COMPILER)
    /// Queue-derived physical grid; logical blocks are traversed with a grid-stride loop.
    std::uint32_t grid_blocks = crb::impl::get_grid_stride_blocks(stream, gpu_kernel::threshold);
    /// Metadata copies retained until the owning queue reports completion.
    crb::DeviceMemory<gpu_kernel::CBCD_Launch> device_launches;
    /// Number of allocated metadata records; storage grows only with an idle lane.
    std::size_t device_launch_capacity = 0;
#endif
    /// By-value launches of the next submission. Only the first num_launches are used.
    std::vector<gpu_kernel::CBCD_Launch> launches;
    std::size_t num_launches = 0;
    std::vector<std::uint64_t> destination_values;
    std::vector<std::size_t> destination_peers;
    // Every appended face aggregate supplies all four fields in AddBatchSlot.
    struct OutgoingFace // NOLINT(cppcoreguidelines-pro-type-member-init)
    {
      /// Index of the launch containing the payload offset.
      std::uint32_t launch;
      /// Index of the uint64 payload-offset entry within that launch.
      std::uint32_t data_entry;
      /// Outgoing-bank range index, also indexing the communicator's destination peers.
      std::size_t destination_peer;
      /// The face as published once the launch completes; its offset is set at submission.
      CBCD_OutgoingFaceRecord record;
    };
    std::vector<OutgoingFace> outgoing_faces;
    bool in_flight = false;
    bool completed = false;

    void Poll()
    {
      completed = in_flight and stream.is_completed();
      if (completed)
        in_flight = false;
    }
  };

  struct AngleSetLaunch
  {
    /// Threads per block along the angle-group (x) and cell (y) axes.
    crb::Dim3 threads_per_block;
    /// Number of blocks spanning the angle-group axis.
    unsigned int num_stride_blocks = 0;
    /// Borrowed FLUDS, held by the angle set through a shared pointer.
    CBCD_FLUDS* fluds = nullptr;
  };

  DiscreteOrdinatesProblem& problem_;
  /// Device grid limit in logical launch blocks. CUDA and HIP use one grid block per launch block.
  std::uint32_t max_launch_blocks_;
  /// Launch functions by metadata size class, selected for the rank's local cell DOF counts.
  std::array<gpu_kernel::CBCD_LaunchFunction, gpu_kernel::CBCD_NUM_LAUNCH_CLASSES>
    launch_functions_{};
  /// True when a boundary has lagged angular fluxes for this groupset.
  bool needs_boundary_upload_ = false;
  std::unique_ptr<CBCD_AsynchronousCommunicator> async_comm_;
  /// Borrowed angle sets, owned by the groupset's angle aggregation.
  std::vector<CBCD_AngleSet*> angle_sets_;
  std::vector<AngleSetLaunch> angle_set_launches_;
  /// Host copy of the sweep arguments indexed by angle set.
  crb::HostVector<gpu_kernel::Arguments<SweepKind::CBC>> host_arguments_;
  /// Device copy of the sweep arguments indexed by angle set.
  crb::DeviceMemory<gpu_kernel::Arguments<SweepKind::CBC>> device_arguments_;
  /// Device table of FLUDS-owned saved-flux pointers indexed by angle set.
  /// Entries remain valid until the chunk is rebuilt and are nullptr when saving is disabled.
  crb::DeviceMemory<double*> device_saved_psi_;
  std::vector<std::unique_ptr<FusedLane>> lanes_;
};

} // namespace opensn
