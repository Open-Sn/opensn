// SPDX-FileCopyrightText: 2026 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep/communicators/lock_free_queues.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep/scheduler/spmd_threadpool.h"
#include <atomic>
#include <cassert>
#include <cstddef>
#include <cstdint>
#include <memory>
#include <vector>

namespace opensn
{

class AngleSet;
class SweepCommunicator;

/**
 * Return the sweep worker that owns an angle set.
 * The ASYNC_FIFO scheduler assigns the angle sets of a groupset to its workers cyclically. The
 * communicator sizes its single-producer queues for the same assignment.
 *
 * \param angle_set_id Angle-set index.
 * \param num_workers Number of sweep workers.
 * \return Index of the owning worker.
 */
constexpr std::size_t
CBCD_AngleSetWorker(std::size_t angle_set_id, std::size_t num_workers)
{
  return angle_set_id % num_workers;
}

/// One outgoing nonlocal face published by a sweep worker.
struct CBCD_OutgoingFaceRecord
{
  /// Index of the angle set that produced the face.
  std::size_t angle_set_id = 0;
  /// Index of the face among the faces that the destination rank receives from this angle set.
  std::uint32_t destination_face_index = 0;
  /// Device angular fluxes of the face in destination-node order, valid until the sweep ends.
  const double* psi_values = nullptr;
  /// Number of angular-flux values of the face.
  std::size_t num_psi_values = 0;
};

/**
 * Outgoing nonlocal faces of a groupset's angle sets, queued by destination and sweep worker.
 * Each (destination, worker) pair has its own single-producer, single-consumer queue. The worker
 * that owns an angle set produces its faces, and the communicator's progress thread consumes
 * them. Every queue holds all of its faces of one sweep, so a producer never waits.
 */
class CBCD_OutgoingQueues
{
public:
  using Queue = LockFreeSPSCSlotQueue<CBCD_OutgoingFaceRecord>;

  /// Outgoing queues toward one destination.
  struct Channel
  {
    /// Rank of the destination in the sweep communicator.
    int peer_rank = 0;
    /// One queue per sweep worker. Queues of workers without faces for this destination are empty.
    std::vector<std::unique_ptr<Queue>> worker_queues;
    /// Workers with at least one outgoing face toward this destination.
    std::vector<std::size_t> active_workers;
  };

  /**
   * Set the destinations and the number of faces each angle set sends to each of them.
   *
   * \param peer_ranks Sweep-communicator ranks of the destinations, indexed by channel.
   * \param channels_by_angle_set For each angle set, its channels indexed by the angle set's
   *        destination index.
   * \param faces_by_angle_set For each angle set, its number of faces per sweep indexed by the
   *        angle set's destination index.
   */
  void Initialize(const std::vector<int>& peer_ranks,
                  std::vector<std::vector<std::size_t>> channels_by_angle_set,
                  std::vector<std::vector<std::size_t>> faces_by_angle_set);

  /**
   * Size the worker queues for a number of sweep workers.
   * The queues are rebuilt only when the number of workers changes.
   *
   * \param num_workers Number of sweep workers.
   */
  void Configure(std::size_t num_workers);

  /**
   * Publish one completed outgoing nonlocal face.
   * Only the worker that owns the angle set may call this during a sweep.
   *
   * \param destination_index Index of the destination in the angle set's destination list.
   * \param worker_id Sweep worker index.
   * \param angle_set_id Angle-set index.
   * \param destination_face_index Index of the face among the destination's faces.
   * \param psi_values Device angular fluxes of the face, in destination-node order.
   * \param num_psi_values Number of angular-flux values.
   */
  void Enqueue(std::size_t destination_index,
               std::size_t worker_id,
               std::size_t angle_set_id,
               std::uint32_t destination_face_index,
               const double* psi_values,
               std::size_t num_psi_values)
  {
    assert(angle_set_id < channels_by_angle_set_.size());
    const auto& channels = channels_by_angle_set_[angle_set_id];
    assert(destination_index < channels.size());
    auto& queue = *channels_[channels[destination_index]].worker_queues[worker_id];
    auto& record = queue.ReserveSlot();
    record.angle_set_id = angle_set_id;
    record.destination_face_index = destination_face_index;
    record.psi_values = psi_values;
    record.num_psi_values = num_psi_values;
    queue.PublishSlot();
  }

  /// Return the channels indexed by destination.
  std::vector<Channel>& GetChannels() { return channels_; }

  /// Return the channels indexed by destination.
  const std::vector<Channel>& GetChannels() const { return channels_; }

  /// Return true when every queue is empty.
  bool AllEmpty() const;

private:
  /// Channels indexed by destination.
  std::vector<Channel> channels_;
  /// For each angle set, its channels indexed by the angle set's destination index.
  std::vector<std::vector<std::size_t>> channels_by_angle_set_;
  /// For each angle set, its number of faces per sweep indexed by its destination index.
  std::vector<std::vector<std::size_t>> faces_by_angle_set_;
  /// Number of sweep workers the queues are sized for (0 before the first Configure).
  std::size_t num_workers_ = 0;
};

/**
 * Interprocess communicator for CBC sweeps on GPU devices.
 * One communicator serves every angle set of a groupset. Sweep workers publish completed outgoing
 * faces to single-producer queues. A progress thread packs them into device packets, exchanges
 * the packets with GPU-aware MPI, and notifies the receiving angle sets through single-consumer
 * mailboxes.
 */
class CBCD_AsynchronousCommunicator
{
public:
  /**
   * Construct the communicator of a groupset's angle sets.
   * The queue capacities and buffer sizes are derived from the angle sets' FLUDS common data.
   * Ranks that exchange angular fluxes require GPU-aware MPI, which a multi-rank problem declares
   * by configuring OpenSn with OPENSN_USE_GPU_AWARE_MPI. This throws when MPI thread support is
   * below MPI_THREAD_SERIALIZED.
   *
   * \param angle_sets Angle sets served by this communicator. Must not be empty.
   * \param sweep_communicator Sweep communication context.
   */
  CBCD_AsynchronousCommunicator(const std::vector<AngleSet*>& angle_sets,
                                const SweepCommunicator& sweep_communicator);

  ~CBCD_AsynchronousCommunicator();

  /// Return the device buffer of received angular fluxes, valid until destruction.
  double* GetIncomingPacketValues() const;

  /// Publish one completed outgoing nonlocal face (see CBCD_OutgoingQueues::Enqueue).
  void EnqueueOutgoing(std::size_t destination_index,
                       std::size_t worker_id,
                       std::size_t angle_set_id,
                       std::uint32_t destination_face_index,
                       const double* psi_values,
                       std::size_t num_psi_values)
  {
    outgoing_.Enqueue(destination_index,
                      worker_id,
                      angle_set_id,
                      destination_face_index,
                      psi_values,
                      num_psi_values);
  }

  /// Return true when the angle set has received incoming faces that are not yet processed.
  bool HasIncoming(std::size_t angle_set_id) const
  {
    return not incoming_dependencies_[angle_set_id]->Empty();
  }

  /**
   * Process the angle set's received incoming faces.
   * Only the worker that owns the angle set may call this during a sweep.
   *
   * \param angle_set_id Angle-set index.
   * \param callback Called with the local ID of the downwind cell of each received face.
   * \return True when at least one face was processed.
   */
  template <typename Callback>
  bool ProcessIncomingDependencies(std::size_t angle_set_id, Callback callback)
  {
    return incoming_dependencies_[angle_set_id]->ProcessReady(callback) > 0;
  }

  /// Mark an angle set complete after all of its outgoing faces have been published.
  void SignalAngleSetComplete(std::size_t angle_set_id);

  /**
   * Start a sweep on the progress thread.
   * The worker queues are rebuilt when the number of workers changes.
   *
   * \param num_workers Number of sweep workers.
   */
  void Start(std::size_t num_workers);

  /// Wait until every angle set is complete and all of the sweep's data has been exchanged.
  void Stop();

private:
  struct PacketExchange;

  /// Progress-thread body: exchange packets until Stop is requested and the sweep is drained.
  void CommThreadLoop();
  /// Return true when every angle set is complete and every outgoing queue is empty.
  bool AllAngleSetsComplete() const;

  const SweepCommunicator& sweep_communicator_;
  /// Outgoing queues of the sweep workers.
  CBCD_OutgoingQueues outgoing_;
  /// Mailboxes of received downwind-cell IDs from the progress thread, indexed by angle set.
  std::vector<std::unique_ptr<LockFreeSPSCSlotQueue<std::uint32_t>>> incoming_dependencies_;
  /// Device buffers, MPI requests, and progress state of the packet exchange.
  std::unique_ptr<PacketExchange> packet_exchange_;
  /// Set by Stop to let the progress thread exit once the sweep is drained.
  std::atomic<bool> stop_requested_{false};
  /// Completion flags indexed by angle set.
  std::vector<std::atomic<bool>> angle_set_complete_;
  /// True between Start and Stop.
  bool sweep_active_ = false;
  /// Pool holding the progress thread, reused across sweeps.
  SPMD_ThreadPool comm_pool_;
};

} // namespace opensn
