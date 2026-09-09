// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep/spds/spds.h"
#include <cstdint>
#include <unordered_map>
#include <utility>

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
   * \param face_neighbor_info Cached neighbor information for every local cell face.
   * \param allow_cycles Whether cycles are allowed in the local and global sweep dependency graphs.
   * \param use_gpus Whether to allocate device memory for GPU acceleration.
   */
  AAH_SPDS(int id,
           const Vector3& omega,
           std::shared_ptr<MeshContinuum> grid,
           const SPDSFaceNeighborInfoVec& face_neighbor_info,
           bool allow_cycles,
           bool use_gpus = false);

  /// Returns the id of this SPDS.
  int GetId() const { return id_; }

  /// Returns true if this SPDS is allowed to remove local and global sweep cycles.
  bool AllowCycles() const { return allow_cycles_; }

  /// Return this location's depth in the global sweep graph for DOG scheduling.
  int GetLocationDepth() const { return location_depth_; }

  /// Builds the Feedback Arc Set (FAS) for the global sweep.
  void BuildGlobalSweepFAS();

  /// Builds the global sweep TDG. May be called only once.
  void BuildGlobalSweepTDG();

  /// Copies the levelized SPLS data on device.
  void CopySPLSDataOnDevice();

  /// Free the memory on GPU.
  void FreeDeviceData();

  /// Get level vector on device
  std::uint32_t* GetDeviceLevelVector(std::size_t level) const
  {
    return device_levelized_spls_ + contiguous_offset_[level];
  }

  /// Returns the global sweep FAS as a vector of edges.
  const std::vector<int>& GetGlobalSweepFAS() const { return global_sweep_fas_; }

  /**
   * Sets the global sweep FAS.
   * \param edges A vector of edges representing the FAS.
   */
  void SetGlobalSweepFAS(std::vector<int> edges) { global_sweep_fas_ = std::move(edges); }

  /// Returns non-zero outgoing edge weights as (destination, weight) pairs.
  std::vector<std::pair<int, double>> ComputeLocalLocationEdgeWeights() const;

  /// Sets sparse edge weights keyed by `dep * comm_size + loc`.
  void SetGlobalEdgeWeights(std::unordered_map<std::int64_t, double> weights)
  {
    global_edge_weights_ = std::move(weights);
  }

  /// Sets the global location-to-location dependencies (result of BatchCommunicateLocationDeps).
  void SetGlobalDependencies(std::vector<std::vector<int>> deps)
  {
    global_dependencies_ = std::move(deps);
  }

  /// Destructor.
  ~AAH_SPDS() override;

private:
  /// Unique identifier for this SPDS.
  int id_;
  /// Flag indicating whether cycles are allowed in the dependency graphs.
  bool allow_cycles_;
  /// Location-to-location global sweep dependencies.
  std::vector<std::vector<int>> global_dependencies_;
  /// This location's depth in the global sweep graph for DOG scheduling.
  int location_depth_ = -1;
  /// Whether the global sweep TDG has been built and its setup data released.
  bool global_sweep_tdg_built_ = false;
  /// Vector of edges representing the FAS used to break cycles in the global sweep graph.
  std::vector<int> global_sweep_fas_;
  /// Sparse global edge weights, keyed by `dep * comm_size + loc`.
  std::unordered_map<std::int64_t, double> global_edge_weights_;
  /// Levelized SPLS structure on GPU (only visible to GPU implementation).
  std::uint32_t* device_levelized_spls_ = nullptr;
  /// Per-level offset into the contiguous level data.
  std::vector<std::uint64_t> contiguous_offset_;
};

} // namespace opensn
