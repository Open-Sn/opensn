#pragma once

#include <utility>
#include "framework/mesh/VolumeMesher/chi_volumemesher.h"
#include "framework/mesh/Cell/cell.h"

namespace chi_mesh
{

/**
 * An extruder mesher taking a flat surface and extruding it.
 */
class VolumeMesherExtruder : public VolumeMesher
{
public:
  enum class TemplateType : int
  {
    UNPARTITIONED_MESH = 2
  };
  struct MeshLayer
  {
    std::string name;
    double height;
    int sub_divisions;
  };

private:
  const TemplateType template_type_;
  std::shared_ptr<const UnpartitionedMesh> template_unpartitioned_mesh_ = nullptr;

  std::vector<MeshLayer> input_layers_;
  std::vector<double> vertex_layers_;
  size_t node_z_index_incr_ = 0;

  uint64_t zmax_bndry_id = 4;
  uint64_t zmin_bndry_id = 5;

public:
  explicit VolumeMesherExtruder(std::shared_ptr<const UnpartitionedMesh> in_unpartitioned_mesh)
    : VolumeMesher(VolumeMesherType::EXTRUDER),
      template_type_(TemplateType::UNPARTITIONED_MESH),
      template_unpartitioned_mesh_(std::move(in_unpartitioned_mesh))
  {
  }

  const std::vector<double>& GetVertexLayers() const { return vertex_layers_; }
  void AddLayer(const MeshLayer& new_layer) { input_layers_.push_back(new_layer); }

  void Execute() override;

private:
  /**
   * Creates actual z-levels for the input layer specification.
   */
  void SetupLayers(int default_layer_count = 1);

  /**
   * Projects a centroid to an extruded equivalent layer.
   */
  Vector3 ProjectCentroidToLevel(const Vector3& centroid, size_t level);

  /**
   * Computes a cell's partition id based on a centroid.
   */
  int GetCellKBAPartitionIDFromCentroid(Vector3& centroid);

  /**
   * Determines if a template cell is in the current partition or a direct neighbor to the current
   * partition.
   */
  bool
  HasLocalScope(const Cell& template_cell, const MeshContinuum& template_continuum, size_t z_level);

  /**
   * Makes an extruded cell from a template cell.
   */
  std::unique_ptr<Cell> MakeExtrudedCell(const Cell& template_cell,
                                         const MeshContinuum& grid,
                                         size_t z_level,
                                         uint64_t cell_global_id,
                                         int partition_id,
                                         size_t num_template_cells);

  /**
   * Creates nodes that are owned locally from the 2D template grid.
   */
  void CreateLocalNodes(MeshContinuum& template_grid, MeshContinuum& grid);

  /**
   * Extrude template cells into polygons.
   */
  void ExtrudeCells(MeshContinuum& template_grid, MeshContinuum& grid);
};

} // namespace chi_mesh
