// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include "framework/object.h"
#include "framework/mesh/unpartitioned_mesh/unpartitioned_mesh.h"
#include "mpicpp-lite/mpicpp-lite.h"

namespace mpi = mpicpp_lite;

namespace opensn
{
class GraphPartitioner;
class MeshContinuum;

/**
 * Mesh generation can be very complicated in parallel. Some mesh formats
 * do not store connectivity information and therefore we have to establish
 * connectivity after the file is read. Most mesh formats are also not
 * partitioned for multiple processors when read/generated. The design of these
 * types of objects need to consider the previous generation of how we did this.
 * We had the concept of a VolumeMesher, which essentially took an
 * `UnpartitionedMesh` and converted it into a partitioned mesh (i.e. the
 * required `MeshContinuum` object). Up to now we've really only had two
 * variants. The `VolumeMesherPredefinedUnpartitioned` and the
 * `VolumeMesherExtruder`. Both of which operated on an `UnpartitionedMesh`
 * object. With this new design we want to unify these concepts to make them
 * more extendable and therefore we split a `MeshGenerator`'s execution into a
 * phase that generates an unpartitioned mesh and a phase that then converts
 * this mesh into real mesh (with both steps customizable). The phase that
 * creates the real mesh can be hooked up to a partitioner that can also be
 * designed to be pluggable.
 */
class MeshGenerator : public Object
{
public:
  /**
   * Final execution step.
   */
  virtual void Execute();

  static InputParameters GetInputParameters();
  explicit MeshGenerator(const InputParameters& params);

  /**
   * Virtual method to generate the unpartitioned mesh for the next step.
   *
   * Default behavior here is to return the input umesh unaltered.
   */
  virtual std::shared_ptr<UnpartitionedMesh>
  GenerateUnpartitionedMesh(std::shared_ptr<UnpartitionedMesh> input_umesh);

  struct VertexListHelper
  {
    virtual const Vertex& at(uint64_t vid) const = 0;
  };
  template <typename T>
  struct STLVertexListHelper : public VertexListHelper
  {
    explicit STLVertexListHelper(const T& list) : list_(list) {}
    const Vertex& at(uint64_t vid) const override { return list_.at(vid); };
    const T& list_;
  };

protected:
  /**
   * Builds a cell-graph and executes the partitioner that assigns cell
   * partition ids based on the supplied number of partitions.
   */
  std::vector<int64_t> PartitionMesh(const UnpartitionedMesh& input_umesh, int num_partitions);

  /**
   * Executes the partitioner and configures the mesh as a real mesh.
   */
  std::shared_ptr<MeshContinuum> SetupMesh(std::shared_ptr<UnpartitionedMesh> input_umesh,
                                           const std::vector<int64_t>& cell_pids);

  /**
   * Broadcasts PIDs to other locations.
   */
  static void
  BroadcastPIDs(std::vector<int64_t>& cell_pids, int root, const mpi::Communicator& communicator);

  /**
   * Determines if a cells needs to be included as a ghost or as a local cell.
   */
  bool CellHasLocalScope(int location_id,
                         const UnpartitionedMesh::LightWeightCell& lwcell,
                         uint64_t cell_global_id,
                         const std::vector<std::set<uint64_t>>& vertex_subscriptions,
                         const std::vector<int64_t>& cell_partition_ids) const;

  /**
   * Converts a light-weight cell to a real cell.
   */
  static std::unique_ptr<Cell> SetupCell(const UnpartitionedMesh::LightWeightCell& raw_cell,
                                         uint64_t global_id,
                                         uint64_t partition_id,
                                         const VertexListHelper& vertices);

  static void SetGridAttributes(MeshContinuum& grid,
                                MeshAttributes new_attribs,
                                std::array<size_t, 3> ortho_cells_per_dimension);

  static void ComputeAndPrintStats(const MeshContinuum& grid);

  const double scale_;
  const bool replicated_;
  std::vector<MeshGenerator*> inputs_;
  GraphPartitioner* partitioner_ = nullptr;
};

} // namespace opensn
