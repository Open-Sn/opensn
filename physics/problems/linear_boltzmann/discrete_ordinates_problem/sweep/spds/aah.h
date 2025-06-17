// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include "physics/problems/linear_boltzmann/discrete_ordinates_problem/sweep/spds/spds.h"

namespace opensn
{

class AAH_SPDS : public SPDS
{
public:
  /**
   * Creates a sweep-plane data structure (SPDS) for the given direction and grid.
   *
   * \param id The unique identifier for this SPDS.
   * \param omega The angular direction for the sweep operation.
   * \param grid The grid on which the sweep is performed.
   * \param allow_cycles Whether cycles are allowed in the local and global swepp dependency graphs.
   */
  AAH_SPDS(int id,
           const Vector3& omega,
           const std::shared_ptr<MeshContinuum> grid,
           bool allow_cycles);

  /// Returns the id of this SPDS.
  int GetId() { return id_; }

  /// Return the levelized global sweep TDG.
  const std::vector<STDG>& GetGlobalSweepPlanes() const { return global_sweep_planes_; }

  /// Builds the Feedback Arc Set (FAS) for the global sweep.
  void BuildGlobalSweepFAS();

  /// Builds the Task Dependency Graph (TDG) for the global sweep.
  void BuildGlobalSweepTDG();

  /// Returns the global sweep FAS as a vector of edges.
  std::vector<int> GetGlobalSweepFAS() { return global_sweep_fas_; }

  /**
   * Sets the global sweep FAS.
   * \param edges A vector of edges representing the FAS.
   */
  void SetGlobalSweepFAS(std::vector<int>& edges) { global_sweep_fas_ = edges; }

private:
  /// Unique identifier for this SPDS.
  int id_;
  /// Flag indicating whether cycles are allowed in the dependency graphs.
  bool allow_cycles_;
  /// Location-to-location global sweep dependencies.
  std::vector<std::vector<int>> global_dependencies_;
  /// Levelized global sweep task dependency graph.
  std::vector<STDG> global_sweep_planes_;
  /// Vector of edges representing the FAS used to break cycles in the global sweep graph.
  std::vector<int> global_sweep_fas_;
};

} // namespace opensn
