// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep/spds/spds.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep/sweep.h"
#include <cstddef>
#include <cstdint>
#include <limits>
#include <vector>

namespace opensn
{

/**
 * Cell-by-cell sweep plane data structure and exact local-face slot metadata.
 *
 * This class stores the local CBC task DAG for one sweep direction, its topological ordering,
 * and the local directed-face metadata consumed by both the host CBC FLUDS and the device
 * CBCD FLUDS. Each outgoing local face is represented by one directed-face task with a
 * producer-cell rank, a consumer-cell rank, and a face-node count.
 *
 * `ComputeMaxNumLocalPsiSlots()` computes the exact minimum safe number of reusable local-face
 * slots and the corresponding static face-to-slot map. The result is shared by the host and
 * device CBC implementations. `UpdateLocalFaceSlotLayout()` converts the face-to-slot mapping
 * plan into a compact slot bank sized by the maximum face-node extent within each slot.
 */
class CBC_SPDS : public SPDS
{
public:
  /// Value returned when a local face does not participate in the requested face-task map.
  static constexpr std::uint32_t INVALID_LOCAL_FACE_TASK_ID =
    std::numeric_limits<std::uint32_t>::max();

  /// Construct the CBC sweep plane data structure for one angular direction.
  CBC_SPDS(const Vector3& omega, const std::shared_ptr<MeshContinuum>& grid, bool allow_cycles);

  /// Return the local CBC task list.
  const std::vector<Task>& GetTaskList() const noexcept;

  /**
   * Compute the exact minimum number of reusable local-face psi slots.
   *
   * The local directed faces define a poset under the safe reuse relation. A chain in this
   * poset is one statically reusable slot. The required slot count is therefore the minimum
   * chain-cover cardinality, equivalently the maximum antichain cardinality by Dilworth's
   * theorem.
   *
   * The planner obtains this value from a maximum cardinality matching on the bipartite
   * split graph of the reuse relation. The extracted chains define the static slot IDs.
   * `UpdateLocalFaceSlotLayout()` then sizes each slot by the maximum face-node extent
   * over the faces assigned to that slot.
   */
  void ComputeMaxNumLocalPsiSlots();

  std::size_t GetMaxNumLocalPsiSlots() const noexcept { return max_num_local_psi_slots_; }

  const std::vector<std::uint32_t>& GetLocalFaceSlotIDs() const noexcept
  {
    return local_face_slot_ids_;
  }

  const std::vector<std::uint32_t>& GetLocalFaceSlotNodeOffsets() const noexcept
  {
    return local_face_slot_node_offsets_;
  }

  const std::vector<std::uint16_t>& GetLocalFaceSlotNodeCounts() const noexcept
  {
    return local_face_slot_node_counts_;
  }

  std::size_t GetTotalLocalFaceSlotNodes() const noexcept { return total_local_face_slot_nodes_; }

  std::size_t GetMaxLocalFaceNodeCount() const noexcept { return max_local_face_node_count_; }

  /// Return the local directed-face task ID for an outgoing local face.
  std::uint32_t GetOutgoingLocalFaceTaskID(std::uint32_t cell_local_id,
                                           unsigned int face_id) const noexcept;

  /// Return the local directed-face task ID for an incoming local face.
  std::uint32_t GetIncomingLocalFaceTaskID(std::uint32_t cell_local_id,
                                           unsigned int face_id) const noexcept;

  ~CBC_SPDS() override = default;

private:
  /// Build the local cell task DAG and its successor-rank adjacency.
  void BuildTaskGraph();

  /// Enumerate local directed faces and map them to producer and consumer cell ranks.
  void BuildLocalFaceTaskGraph();

  /// Topological ordering of local cell IDs: topo_order_[rank] = cell_local_id.
  std::vector<std::uint32_t> topo_order_;
  /// Topological rank keyed by local cell ID.
  std::vector<std::uint32_t> topo_rank_by_cell_local_id_;
  /// Per-cell task descriptors with successor adjacency lists.
  std::vector<Task> task_list_;
  /// Offsets into the flat successor-rank array indexed by topological task rank.
  std::vector<std::uint32_t> task_successor_rank_offsets_;
  /// Flat successor topological ranks grouped by producer task rank.
  std::vector<std::uint32_t> task_successor_ranks_;
  /// Flat face-table offsets indexed by cell local IDs.
  std::vector<std::uint32_t> cell_face_offsets_;
  /// Flat outgoing local-face task IDs indexed by face storage index.
  std::vector<std::uint32_t> outgoing_local_face_task_ids_;
  /// Flat incoming local-face task IDs indexed by face storage index.
  std::vector<std::uint32_t> incoming_local_face_task_ids_;
  /// Face-rank offsets grouped by producer-cell topological rank.
  std::vector<std::uint32_t> producer_cell_face_offsets_;
  /// Producer-cell topological rank for each local directed face.
  std::vector<std::uint32_t> local_face_producer_ranks_;
  /// Consumer-cell topological rank for each local directed face.
  std::vector<std::uint32_t> local_face_consumer_ranks_;
  /// Number of nodes for each local directed face task.
  std::vector<std::uint16_t> local_face_node_counts_;
  /// Static slot assignment: local_face_slot_ids_[face_task_id] = slot_id.
  std::vector<std::uint32_t> local_face_slot_ids_;
  /// Slot-local node extents: local_face_slot_node_counts_[slot_id] = max nodes in that slot.
  std::vector<std::uint16_t> local_face_slot_node_counts_;
  /// Prefix offsets into the compact local-face slot bank.
  std::vector<std::uint32_t> local_face_slot_node_offsets_;
  /// Minimum number of local-face angular flux storage slots.
  std::size_t max_num_local_psi_slots_ = 0;
  /// Total number of local-face nodes in the compact slot bank.
  std::size_t total_local_face_slot_nodes_ = 0;
  /// Maximum number of nodes across all local directed faces.
  std::size_t max_local_face_node_count_ = 0;

  /**
   * Recompute slot-local node extents and prefix offsets from the current slot assignment.
   *
   * Each slot is sized to the maximum face-node extent of the local directed faces assigned
   * to that slot. This preserves the exact slot count while avoiding one global slot extent.
   */
  void UpdateLocalFaceSlotLayout();
};

} // namespace opensn
