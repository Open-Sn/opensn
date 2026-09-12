// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep/spds/spds.h"
#include <cstdint>
#include <stdexcept>
#include <utility>
#include <vector>

namespace opensn
{

class AAH_SPDS : public SPDS
{
public:
  struct GlobalSweepMetadata
  {
    std::vector<int> location_depths;
    std::vector<std::vector<int>> delayed_dependencies;
    std::vector<std::vector<int>> delayed_successors;
  };

  struct GlobalSweepEdge
  {
    int dependency;
    int location;
    double weight;
  };

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

  /// Builds the global sweep metadata on the graph-owning rank.
  GlobalSweepMetadata BuildGlobalSweepMetadata();

  /// Installs the metadata needed by this rank.
  void SetGlobalSweepMetadata(int location_depth,
                              std::vector<int> delayed_dependencies,
                              std::vector<int> delayed_successors);

  /// Copies the levelized SPLS data on device.
  void CopySPLSDataOnDevice();

  /// Free the memory on GPU.
  void FreeDeviceData();

  /// Get level vector on device
  std::uint32_t* GetDeviceLevelVector(std::size_t level) const
  {
    return device_levelized_spls_ + contiguous_offset_[level];
  }

  /// Returns non-zero outgoing edge weights as (destination, weight) pairs.
  std::vector<std::pair<int, double>> ComputeLocalLocationEdgeWeights() const;

  /// Sets the sparse global sweep graph edges.
  void SetGlobalEdges(std::vector<GlobalSweepEdge> edges)
  {
    if (global_edges_initialized_ or global_sweep_metadata_built_)
      throw std::logic_error("AAH_SPDS: Global sweep edges have already been set.");
    global_edges_ = std::move(edges);
    global_edges_initialized_ = true;
  }

  /// Destructor.
  ~AAH_SPDS() override;

private:
  /// Unique identifier for this SPDS.
  int id_;
  /// Flag indicating whether cycles are allowed in the dependency graphs.
  bool allow_cycles_;
  /// This location's depth in the global sweep graph for DOG scheduling.
  int location_depth_ = -1;
  bool global_edges_initialized_ = false;
  bool global_sweep_metadata_built_ = false;
  bool global_sweep_metadata_set_ = false;
  std::vector<GlobalSweepEdge> global_edges_;
  /// Levelized SPLS structure on GPU (only visible to GPU implementation).
  std::uint32_t* device_levelized_spls_ = nullptr;
  /// Per-level offset into the contiguous level data.
  std::vector<std::uint64_t> contiguous_offset_;
};

} // namespace opensn
