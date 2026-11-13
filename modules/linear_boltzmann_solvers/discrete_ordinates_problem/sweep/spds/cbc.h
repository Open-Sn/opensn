// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep/spds/spds.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep/spds/cbc_host_task_graph.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep/sweep.h"
#include <cstddef>
#include <cstdint>
#include <limits>
#include <set>
#include <string>
#include <span>
#include <unordered_map>
#include <utility>
#include <vector>

namespace opensn
{

/// CBC sweep-plane data structure.
class CBC_SPDS : public SPDS
{
public:
  /// Backend-owned task representation; five-argument construction retains the legacy task API.
  enum class TaskStorage
  {
    DEVICE_TASKS,
    HOST_CSR
  };

  /// Weighted directed edge in the interpartition CBC sweep graph.
  struct LocationEdgeWeight
  {
    /// Upstream MPI rank.
    int upstream_location = 0;
    /// Downstream MPI rank.
    int downstream_location = 0;
    /// Accumulated sweep-graph edge weight.
    double weight = 0.0;
  };

  /// Directed dependency through one local face, from the upwind to the downwind cell.
  struct LocalFaceEdge
  {
    /// Local ID of the upwind cell.
    std::uint32_t upwind = 0;
    /// Local ID of the downwind cell.
    std::uint32_t downwind = 0;
    /// Face weight, the face-node count.
    double weight = 0.0;
  };

  /**
   * Add the local cell graph and remove a weighted feedback arc set (FAS) from it.
   *
   * Faces between the same ordered pair of cells become one edge whose weight is the sum of their
   * weights, added in the order of the first such face, so the FAS lists every removed cell pair
   * once. Lagging a removed pair lags all faces between the two cells. The shared algorithm is
   * exact for small SCCs and greedy above MAX_EXACT_FAS_SCC_SIZE.
   *
   * \param graph Graph with one vertex per local cell and no edges. On return it holds the
   *        acyclic cell graph.
   * \param face_edges One entry per outgoing local face. The upwind and downwind cells differ.
   * \return The removed (upwind, downwind) cell pairs.
   */
  static std::vector<std::pair<std::uint32_t, std::uint32_t>>
  RemoveLocalCellCycles(Graph& graph, const std::vector<LocalFaceEdge>& face_edges);

  /// Slot node offset of a cell face that is not a non-delayed local face.
  static constexpr std::uint32_t INVALID_SLOT_NODE_OFFSET =
    std::numeric_limits<std::uint32_t>::max();

  /**
   * Construct a cell-by-cell sweep-plane data structure for the given direction and grid.
   *
   * \param id Process-independent CBC SPDS ordinal.
   * \param omega The angular direction vector.
   * \param grid Reference to the grid.
   * \param face_neighbor_info Cached neighbor information for every local cell face.
   * \param allow_cycles Whether cycles are allowed in the local sweep dependency graph.
   * \param task_storage Backend task representation. Runtime construction derives this from
   *        use_gpus; the default preserves the existing five-argument C++ constructor contract.
   */
  CBC_SPDS(int id,
           const Vector3& omega,
           const std::shared_ptr<Mesh>& grid,
           const SPDSFaceNeighborInfoVec& face_neighbor_info,
           bool allow_cycles,
           TaskStorage task_storage = TaskStorage::DEVICE_TASKS);

  /// Return the process-independent CBC SPDS ordinal.
  int GetId() const noexcept { return id_; }

  /// Return whether cyclic dependencies may be lagged.
  bool AllowCycles() const noexcept { return allow_cycles_; }

  /// Collectively return a global cell-cycle witness for this no-lag ordering, or an empty string.
  /// All ranks participate, including ranks without cells. No sweep state is modified.
  std::string CheckGlobalCellCycles() const;

  /// Return legacy/device tasks; throws std::logic_error for HOST_CSR storage.
  const std::vector<Task>& GetTaskList() const;

  /**
   * Prepare the host CSR during angle-set construction, before any sweep advances.
   * Legacy/device SPDS callers retain their task API and get one cached CSR adapter
   * shared by host angle sets. HOST_CSR runtime construction already owns the CSR
   * and never materializes legacy tasks. Setup only: do not call concurrently with
   * task rebuilds or sweeps.
   *
   * \return Borrowed immutable CSR until the next task rebuild or SPDS destruction.
   */
  const detail::CBCHostTaskGraph& PrepareHostTaskGraph() const;

  /// Borrow immutable host CSR until task rebuild or SPDS destruction.
  /// Throws std::logic_error if legacy storage has not been prepared for host use.
  /// Reset angle sets before a rebuild.
  const detail::CBCHostTaskGraph& GetHostTaskGraph() const;

  /// Return the initial unsatisfied dependency count indexed by local task ID.
  const std::vector<unsigned int>& GetInitialTaskDependencyCounts() const noexcept
  {
    return initial_task_dependency_counts_;
  }

  /// Return the initially ready local task stack in local-cell order.
  const std::vector<std::uint32_t>& GetInitialReadyTasks() const noexcept
  {
    return initial_ready_tasks_;
  }

  /// Return flattened rank pairs removed from the interpartition sweep graph.
  const std::vector<int>& GetGlobalSweepFAS() const { return global_sweep_fas_; }

  /// Set flattened rank pairs removed from the interpartition sweep graph.
  void SetGlobalSweepFAS(std::vector<int> edges);

  /// Build the global feedback arc set from the interpartition sweep graph.
  void BuildGlobalSweepFAS();

  /// Apply the global feedback arc set to location dependencies and rebuild tasks.
  void ApplyGlobalSweepFAS();

  /**
   * Compute sparse local contributions to interpartition edge weights.
   *
   * Each outgoing rank edge is weighted by its total face-node count, which is proportional to the
   * lagged nonlocal angular-flux storage and copy work for this SPDS.
   *
   * \return The local contributions to the global edge weights.
   */
  std::vector<LocationEdgeWeight> ComputeLocalLocationEdgeWeights() const;

  /// Store sparse global edge weights for global feedback arc set construction.
  void SetGlobalEdgeWeights(std::span<const LocationEdgeWeight> edge_weights);

  /// Store the incoming interpartition dependencies for every MPI rank.
  void SetGlobalDependencies(std::vector<std::vector<int>> global_dependencies);

  /// Return whether a local upwind-to-downwind cell dependency is delayed.
  bool IsDelayedLocalDependency(std::uint32_t upwind_local_id,
                                std::uint32_t downwind_local_id) const noexcept;

  /**
   * Return the total number of nodes of the local-face slot bank: the sum over the slots of the
   * node count of each slot's largest face.
   */
  std::size_t GetTotalLocalFaceSlotNodes() const noexcept { return total_local_face_slot_nodes_; }

  /**
   * Return the cell-face prefix offsets shared with host CBC FLUDS common data.
   * The vector object lives as long as this SPDS. The mesh and slot plan must remain unchanged
   * while common data or FLUDS banks constructed from them are in use.
   */
  const std::vector<std::uint32_t>& GetCellFaceOffsets() const noexcept
  {
    return cell_face_offsets_;
  }

  /**
   * Return the local slot-node offsets indexed by cell face, including INVALID_SLOT_NODE_OFFSET.
   * Ownership and lifetime are the same as GetCellFaceOffsets().
   */
  const std::vector<std::uint32_t>& GetLocalFaceSlotNodeOffsets() const noexcept
  {
    return local_face_slot_node_offsets_;
  }

  /**
   * Return the node offset, in the local-face slot bank, of the slot of a non-delayed local face,
   * for both of its cell faces (the producer's outgoing and the consumer's incoming one), or
   * INVALID_SLOT_NODE_OFFSET for any other cell face (boundary, nonlocal, delayed or parallel).
   * The slot holds the face's nodes in the producer's face-node order.
   */
  std::uint32_t GetLocalFaceSlotNodeOffset(std::uint32_t cell_local_id,
                                           unsigned int face_id) const noexcept
  {
    return local_face_slot_node_offsets_[cell_face_offsets_[cell_local_id] + face_id];
  }

protected:
  /// Build local sweep tasks from current local and delayed dependencies.
  void BuildTaskList();

  /// Build checked host CSR directly, without allocating per-cell successor vectors.
  void BuildHostTaskGraph();

  /// Cache a CSR adapter for legacy host callers, preserving their existing task list.
  void BuildLegacyHostTaskGraph() const;

  /**
   * Assign slots of the local-face slot bank to the non-delayed local faces, with the minimum
   * number of slots (detail::BuildMinimumLocalFaceSlotPlan), and set the slot node offsets of the
   * cell faces and the total number of nodes of the bank.
   */
  void BuildLocalFaceSlotPlan(const SPDSFaceNeighborInfoVec& face_neighbor_info);

  /**
   * Number the non-delayed local faces by the topological rank of their producer task and then by
   * the producer's face index, which is the order of the task successors (BuildTaskList), and
   * write each face's ID at both of its cell faces in local_face_slot_node_offsets_.
   *
   * \param face_neighbor_info Cached neighbor information for every local cell face.
   * \param producer_face_offsets Output CSR offsets of the faces indexed by producer task rank.
   * \param consumer_ranks Output consumer task rank indexed by face ID.
   * \param node_counts Output node count indexed by face ID. It is the face's vertex count, which
   *        equals its number of face nodes for the PWLD discretization that the sweeps use.
   */
  void BuildLocalFaceReuseData(const SPDSFaceNeighborInfoVec& face_neighbor_info,
                               std::vector<std::uint32_t>& producer_face_offsets,
                               std::vector<std::uint32_t>& consumer_ranks,
                               std::vector<std::uint16_t>& node_counts);

  /// Process-independent ordinal used to order CBC collectives.
  int id_ = 0;
  /// Whether cyclic dependencies may be broken by lagging fluxes.
  bool allow_cycles_ = false;
  /// Backend selection, fixed for the SPDS lifetime.
  TaskStorage task_storage_;
  /// Legacy/device tasks, unallocated on the host CSR path.
  std::vector<Task> task_list_;
  /// Host CSR arrays; a device-only SPDS never allocates them.
  mutable detail::CBCHostTaskGraph host_task_graph_;
  /// Contiguous initial dependency counts used to reset every angle-set sweep.
  std::vector<unsigned int> initial_task_dependency_counts_;
  /// Initially ready tasks in local-cell order.
  std::vector<std::uint32_t> initial_ready_tasks_;
  /// Incoming interpartition dependencies for each MPI rank.
  std::vector<std::vector<int>> global_dependencies_;
  /// Flattened pairs of rank edges removed from the global sweep graph.
  std::vector<int> global_sweep_fas_;
  /// Sparse transport weights keyed by directed interpartition edge.
  std::unordered_map<std::uint64_t, double> global_edge_weights_;
  /// Delayed local upwind-to-downwind cell dependencies.
  std::set<std::uint64_t> delayed_local_dependency_set_;
  /// MPI-rank-indexed flags for delayed incoming location dependencies.
  std::vector<unsigned char> delayed_location_dependency_flags_;

  /// Prefix offsets into cell-face-indexed arrays.
  std::vector<std::uint32_t> cell_face_offsets_;
  /// Slot node offset indexed by cell face (see GetLocalFaceSlotNodeOffset).
  std::vector<std::uint32_t> local_face_slot_node_offsets_;
  /// Total number of nodes of the local-face slot bank.
  std::size_t total_local_face_slot_nodes_ = 0;
};

} // namespace opensn
