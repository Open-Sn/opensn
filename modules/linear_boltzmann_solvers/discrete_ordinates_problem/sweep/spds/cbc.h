// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep/spds/spds.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep/sweep.h"
#include <cstdint>
#include <set>
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

  /**
   * Construct a cell-by-cell sweep-plane data structure for the given direction and grid.
   *
   * \param id Process-independent CBC SPDS ordinal.
   * \param omega The angular direction vector.
   * \param grid Reference to the grid.
   * \param face_neighbor_info Cached neighbor information for every local cell face.
   * \param allow_cycles Whether cycles are allowed in the local sweep dependency graph.
   */
  CBC_SPDS(int id,
           const Vector3& omega,
           const std::shared_ptr<MeshContinuum>& grid,
           const SPDSFaceNeighborInfoVec& face_neighbor_info,
           bool allow_cycles);

  /// Return the process-independent CBC SPDS ordinal.
  int GetId() const noexcept { return id_; }

  /// Return whether cyclic dependencies may be lagged.
  bool AllowCycles() const noexcept { return allow_cycles_; }

  /**
   * Remove a deterministic feedback arc set from each strongly connected
   * component (SCC).
   *
   * SCCs with 16 or fewer vertices use an exact minimum-weight dynamic program to
   * find the minimum-weight FAS. SCCs with more than 16 vertices use a deterministic
   * weighted-ordering heuristic to find an approximate minimum-weight FAS.
   */
  static std::vector<std::pair<Vertex, Vertex>> RemoveCyclicDependencies(Graph& graph);

  /// Return the local cell task list.
  const std::vector<Task>& GetTaskList() const;

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
  void SetGlobalSweepFAS(std::vector<int> edges) { global_sweep_fas_ = std::move(edges); }

  /// Build the global feedback arc set from the interpartition sweep graph.
  void BuildGlobalSweepFAS();

  /// Apply the global feedback arc set to location dependencies and rebuild tasks.
  void ApplyGlobalSweepFAS();

  /**
   * Compute sparse local contributions to interpartition edge weights.
   *
   * Each outgoing rank edge is weighted by its total face-node count, which is proportional to the
   * lagged nonlocal angular-flux storage and copy work for this SPDS.
   */
  std::vector<LocationEdgeWeight> ComputeLocalLocationEdgeWeights() const;

  /// Store sparse global edge weights for global feedback arc set construction.
  void SetGlobalEdgeWeights(std::span<const LocationEdgeWeight> edge_weights);

  /// Store the incoming interpartition dependencies for every MPI rank.
  void SetGlobalDependencies(std::vector<std::vector<int>> global_dependencies);

  /// Return whether a local upwind-to-downwind cell dependency is delayed.
  bool IsDelayedLocalDependency(std::uint32_t upwind_local_id,
                                std::uint32_t downwind_local_id) const noexcept;

protected:
  /// Approximate a minimum-weight FAS using a deterministic weighted vertex ordering.
  static std::vector<std::pair<Vertex, Vertex>>
  FindApproxMinimumFAS(const Graph& graph, const std::vector<Vertex>& component);

  /// Find an exact minimum-weight FAS by dynamic programming over vertex subsets.
  static std::vector<std::pair<Vertex, Vertex>>
  FindExactMinimumFAS(const Graph& graph, const std::vector<Vertex>& component);

  /// Build local sweep tasks from current local and delayed dependencies.
  void BuildTaskList();

  /// Process-independent ordinal used to order CBC collectives.
  int id_ = 0;
  /// Whether cyclic dependencies may be broken by lagging fluxes.
  bool allow_cycles_ = false;
  /// Cell-by-cell task list.
  std::vector<Task> task_list_;
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
};

} // namespace opensn
