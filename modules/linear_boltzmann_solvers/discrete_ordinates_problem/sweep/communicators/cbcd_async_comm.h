// SPDX-FileCopyrightText: 2026 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep/communicators/lock_free_spsc_queue.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep/scheduler/spmd_threadpool.h"
#include "mpicpp-lite/mpicpp-lite.h"
#include <algorithm>
#include <atomic>
#include <cstddef>
#include <cstdint>
#include <memory>
#include <span>
#include <vector>

namespace mpi = mpicpp_lite;

namespace opensn
{

/// One MPI progress thread per groupset.
inline constexpr std::size_t CBCD_NUM_PROGRESS_THREADS = 1;

class CBCD_FLUDS;
class CBCD_NonLocalBank;
class SweepCommunicator;
struct CBCD_IncomingNonlocalFace;

/// Map an angle set to its sweep worker. The scheduler and communicator use the same assignment.
/// Requires a positive worker count.
constexpr std::size_t
CBCD_AngleSetWorker(std::size_t angle_set_id, std::size_t num_workers)
{
  return angle_set_id % num_workers;
}

/// Reserve one OpenSn thread for MPI progress, with at least one sweep worker and at most
/// one per angle set. Requires at least one angle set.
constexpr std::size_t
CBCD_NumSweepWorkers(std::size_t num_angle_sets, std::size_t num_threads)
{
  return std::min(num_angle_sets,
                  std::max<std::size_t>(num_threads, CBCD_NUM_PROGRESS_THREADS + 1) -
                    CBCD_NUM_PROGRESS_THREADS);
}

/// One outgoing nonlocal face published by a sweep worker.
struct CBCD_OutgoingFaceRecord
{
  std::uint32_t angle_set_id = 0;
  /// Index of the face among the faces that the destination rank receives from this angle set.
  std::uint32_t destination_face_index = 0;
  /// Offset of the face's values in the outgoing bank, in destination-node order.
  std::size_t offset = 0;
  std::size_t num_values = 0;
};

/// Device CBC communicator shared by a groupset's angle sets.
///
/// The progress thread sends contiguous runs directly from the outgoing bank and receives into
/// the incoming bank. Each run uses two tags. The control message has a uint64 payload length
/// in doubles, followed by one uint64 per face (angle-set index in bits 32–63 and receiver's
/// face index in bits 0–31). A payload contains the face values in control order. Same-source,
/// same-tag ordering pairs controls and payloads. Faces for a peer accumulate while that peer
/// has message pairs in flight. Bank offsets count doubles and control offsets count uint64s.
///
/// The sweep communicator and FLUDS must outlive this communicator, which shares ownership
/// of the nonlocal banks. The progress thread handles MPI traffic. Sweep workers may call MPI
/// to abort on failure. Initialize MPI with a compatible thread level and keep it active
/// through destruction.
// Deliberate padding starts the progress-thread state on its own cache line; comm_pool_ stays
// last so its thread is joined before that state is destroyed.
class CBCD_AsynchronousCommunicator // NOLINT(clang-analyzer-optin.performance.Padding)
{
public:
  /**
   * Construct from all angle sets' FLUDS in angle-set order, with shared nonlocal banks.
   *
   * Use the scheduler's worker count from CBCD_NumSweepWorkers.
   */
  CBCD_AsynchronousCommunicator(const std::vector<CBCD_FLUDS*>& fluds_by_angle_set,
                                int groupset_id,
                                const SweepCommunicator& sweep_communicator,
                                std::size_t num_workers);

  CBCD_AsynchronousCommunicator(const CBCD_AsynchronousCommunicator&) = delete;
  CBCD_AsynchronousCommunicator& operator=(const CBCD_AsynchronousCommunicator&) = delete;

  ~CBCD_AsynchronousCommunicator();

  /// Destination peers use the same indices as the outgoing-bank ranges.
  std::size_t GetNumDestinationPeers() const { return destination_peers_.size(); }

  /// Publish a nonlocal face from its owning sweep worker after the producing kernel completes.
  /// Each worker/peer queue can hold all its faces for the sweep.
  void PublishFace(std::size_t destination_peer,
                   std::size_t worker_id,
                   const CBCD_OutgoingFaceRecord& face)
  {
    destination_peers_[destination_peer].worker_queues[worker_id]->Push(face);
  }

  /// Take newly received downwind-cell local IDs on the owning worker. Queue acquisition makes
  /// the corresponding offsets visible after the payload receives complete. The borrowed host
  /// span remains valid until Start resets the queues or the communicator is destroyed.
  std::span<const std::uint32_t> TakeIncomingDependencies(std::size_t angle_set_id)
  {
    return incoming_cell_queues_[angle_set_id]->Take();
  }

  /**
   * Reset the banks and start progress. Call serially with all workers idle.
   * \throw std::logic_error Start is called again before Stop returns.
   */
  void Start();

  /**
   * Wait for all sends and expected receives to drain. Call serially after all workers
   * finish publishing faces.
   */
  void Stop();

private:
  using FaceQueue = LockFreeSPSCQueue<CBCD_OutgoingFaceRecord>;

  struct DestinationPeer
  {
    int peer_rank = 0;
    /// One queue per sweep worker (nullptr for workers without faces for this destination).
    std::vector<std::unique_ptr<FaceQueue>> worker_queues;
    std::vector<std::size_t> active_workers;
  };

  struct MessagePairSend
  {
    /// Index in destination_peers_, with the same numbering as outgoing-bank ranges.
    std::size_t destination_peer = 0;
    mpi::Request control_request;
    mpi::Request payload_request;
  };

  struct SourcePeer
  {
    int peer_rank = 0;
    std::size_t num_faces = 0;
    /// Number of faces not yet announced by a control message in the current sweep.
    std::size_t remaining_faces = 0;
    /// Range of the source's control entries in receive_control_.
    std::size_t control_begin = 0;
    std::size_t control_end = 0;
    /// Offset of the next control message in receive_control_.
    std::size_t control_offset = 0;
  };

  struct PayloadReceive
  {
    mpi::Request request;
    /// Offset of the message pair's first face entry in receive_control_.
    std::size_t first_face_entry = 0;
    std::size_t num_faces = 0;
  };

  void CommThreadLoop();

  bool Flush();

  bool PollSends();

  bool PollControls();

  bool PollReceives();

  void PostControl(std::size_t source_id);

  /// Record each face's offset, then reuse its control entry for the angle-set and cell IDs needed
  /// when the payload completes. Each face is received once per sweep.
  void WriteFaceOffsets(std::size_t source_id,
                        std::size_t first_face_entry,
                        std::size_t num_faces,
                        std::size_t payload_offset,
                        std::size_t num_values);

  bool OutgoingQueuesEmpty() const;

  /// Borrowed communicator context, immutable while the sweep workers run.
  const SweepCommunicator& sweep_communicator_;
  std::vector<DestinationPeer> destination_peers_;
  /// Mailboxes of received downwind-cell IDs from the progress thread, indexed by angle set.
  std::vector<std::unique_ptr<LockFreeSPSCQueue<std::uint32_t>>> incoming_cell_queues_;

  // Read by the progress thread; set at construction.
  mpi::Communicator comm_;
  int control_tag_ = 0;
  int payload_tag_ = 0;
  std::shared_ptr<CBCD_NonLocalBank> outgoing_bank_;
  std::shared_ptr<CBCD_NonLocalBank> incoming_bank_;
  /// FLUDS indexed by angle set.
  std::vector<CBCD_FLUDS*> fluds_by_angle_set_;
  /// For each angle set, its incoming nonlocal faces indexed by source peer.
  std::vector<std::vector<std::span<const CBCD_IncomingNonlocalFace>>> incoming_faces_;
  /// Source peers, in increasing partition order.
  std::vector<SourcePeer> sources_;
  /// Control entries of the sweep's sent and received message pairs.
  std::vector<std::uint64_t> send_control_;
  std::vector<std::uint64_t> receive_control_;
  std::vector<std::size_t> in_flight_;
  std::vector<mpi::Request> control_requests_;

  // Sweep state written only by the progress thread; sending_ starts on a separate cache line.
  alignas(HardwareInterferenceSize) std::vector<MessagePairSend> sending_;
  std::vector<CBCD_OutgoingFaceRecord> ready_;
  std::vector<PayloadReceive> receiving_;
  std::vector<int> control_indices_;
  std::vector<mpi::Status> control_statuses_;
  /// Offset of the next unused uint64 entry in send_control_.
  std::size_t send_control_offset_ = 0;
  std::size_t active_controls_ = 0;
  std::atomic<bool> stop_requested_{false};
  bool started_ = false;
  /// Created on first Start and reused across sweeps. Declared last so destruction joins the
  /// thread before destroying its state.
  SPMD_ThreadPool comm_pool_;
};

} // namespace opensn
