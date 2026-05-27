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
   * Constructs a cell-by-cell sweep-plane data structure (SPDS) with the given direction and grid.
   *
   * \param id Globally unique SPDS identifier.
   * \param omega The angular direction vector.
   * \param grid Reference to the grid.
   * \param allow_cycles Whether cycles are allowed in the local sweep dependency graph.
   */
  CBC_SPDS(int id,
           const Vector3& omega,
           const std::shared_ptr<MeshContinuum>& grid,
           bool allow_cycles);

  int GetId() const noexcept { return id_; }

  /// Return the local cell task list.
  const std::vector<Task>& GetTaskList() const;

  /// Return flattened rank pairs removed from the interpartition sweep graph.
  std::vector<int> GetGlobalSweepFAS() const { return global_sweep_fas_; }

  /// Set flattened rank pairs removed from the interpartition sweep graph.
  void SetGlobalSweepFAS(std::vector<int>& edges) { global_sweep_fas_ = edges; }

  /// Build the global feedback arc set from the interpartition sweep graph.
  void BuildGlobalSweepFAS();

  /// Apply the global feedback arc set to location dependencies and rebuild tasks.
  void ApplyGlobalSweepFAS();

  /// Compute sparse edge weights from this rank to downstream ranks.
  std::vector<LocationEdgeWeight> ComputeLocalLocationEdgeWeights() const;

  /// Store sparse global edge weights for global feedback arc set construction.
  void SetGlobalEdgeWeights(std::span<const LocationEdgeWeight> edge_weights);

  /// Return whether a local upwind-to-downwind cell dependency is delayed.
  bool IsDelayedLocalDependency(std::uint32_t upwind_local_id,
                                std::uint32_t downwind_local_id) const noexcept;

protected:
  /// Build local sweep tasks from current local and delayed dependencies.
  void BuildTaskList();

  /// Globally unique CBC SPDS identifier.
  int id_ = 0;
  /// Whether cyclic dependencies may be broken by lagging fluxes.
  bool allow_cycles_ = false;
  /// Cell-by-cell task list.
  std::vector<Task> task_list_;
  /// Incoming interpartition dependencies for each MPI rank.
  std::vector<std::vector<int>> global_dependencies_;
  /// Flattened pairs of rank edges removed from the global sweep graph.
  std::vector<int> global_sweep_fas_;
  /// Sparse transport weights keyed by directed interpartition edge.
  std::unordered_map<std::uint64_t, double> global_edge_weights_;
  /// Delayed local upwind-to-downwind cell dependencies.
  std::set<std::uint64_t> delayed_local_dependency_set_;
};

} // namespace opensn
