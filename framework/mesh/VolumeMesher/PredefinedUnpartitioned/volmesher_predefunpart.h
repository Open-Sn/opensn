#pragma once

#include "framework/mesh/VolumeMesher/chi_volumemesher.h"
#include "framework/mesh/UnpartitionedMesh/unpartitioned_mesh.h"

namespace chi_mesh
{

/**
 * This volume mesher merely applies a partitioning of an unpartitioned mesh.
 */
class VolumeMesherPredefinedUnpartitioned : public VolumeMesher
{
private:
  std::shared_ptr<const UnpartitionedMesh> umesh_ptr_ = nullptr;

public:
  explicit VolumeMesherPredefinedUnpartitioned(std::shared_ptr<const UnpartitionedMesh> in_umesh)
    : VolumeMesher(VolumeMesherType::UNPARTITIONED), umesh_ptr_(std::move(in_umesh))
  {
  }

  void Execute() override;

  /**
   * Determines if a chi_mesh::UnpartitionedMesh::LightWeightCell is a neighbor to the current
   * partition for ParMETIS-style partitioning.
   * This method loops over the faces of the lightweight cell and determines the partition-id of
   * each the neighbors. If the neighbor has a partition id equal to that of the current process
   * then it means this reference cell is a neighbor.
   */
  static bool CellHasLocalScope(const UnpartitionedMesh::LightWeightCell& lwcell,
                                uint64_t cell_global_id,
                                const std::vector<std::set<uint64_t>>& vertex_subscriptions,
                                const std::vector<int64_t>& cell_partition_ids);

  /**
   * Applies KBA-style partitioning to the mesh.
   */
  static std::vector<int64_t> KBA(const UnpartitionedMesh& umesh);

  /**
   * Applies KBA-style partitioning to the mesh.
   */
  static std::vector<int64_t> PARMETIS(const UnpartitionedMesh& umesh);

  /**
   * Adds a cell to the grid from a light-weight cell.
   */
  static std::unique_ptr<Cell> MakeCell(const UnpartitionedMesh::LightWeightCell& raw_cell,
                                        uint64_t global_id,
                                        uint64_t partition_id,
                                        const std::vector<Vector3>& vertices);
};

} // namespace chi_mesh
