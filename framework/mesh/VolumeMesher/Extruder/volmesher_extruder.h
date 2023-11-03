#pragma once

#include <utility>

#include "mesh/VolumeMesher/chi_volumemesher.h"
#include "mesh/Cell/cell.h"

//###################################################################
/**An extruder mesher taking a flat surface and extruding it.*/
class chi_mesh::VolumeMesherExtruder : public chi_mesh::VolumeMesher
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
  explicit VolumeMesherExtruder(
    std::shared_ptr<const chi_mesh::UnpartitionedMesh> in_unpartitioned_mesh)
    : VolumeMesher(VolumeMesherType::EXTRUDER),
      template_type_(TemplateType::UNPARTITIONED_MESH),
      template_unpartitioned_mesh_(std::move(in_unpartitioned_mesh))
  {
  }

  const std::vector<double>& GetVertexLayers() const { return vertex_layers_; }
  void AddLayer(const MeshLayer& new_layer) { input_layers_.push_back(new_layer); }

  void Execute() override;

private:
  // utils
  void SetupLayers(int default_layer_count = 1);

  chi_mesh::Vector3 ProjectCentroidToLevel(const chi_mesh::Vector3& centroid, size_t level);
  int GetCellKBAPartitionIDFromCentroid(chi_mesh::Vector3& centroid);

  bool HasLocalScope(const chi_mesh::Cell& template_cell,
                     const chi_mesh::MeshContinuum& template_continuum,
                     size_t z_level);

  std::unique_ptr<chi_mesh::Cell> MakeExtrudedCell(const chi_mesh::Cell& template_cell,
                                                   const chi_mesh::MeshContinuum& grid,
                                                   size_t z_level,
                                                   uint64_t cell_global_id,
                                                   int partition_id,
                                                   size_t num_template_cells);

  void CreateLocalNodes(chi_mesh::MeshContinuum& template_grid, chi_mesh::MeshContinuum& grid);

  void ExtrudeCells(chi_mesh::MeshContinuum& template_grid, chi_mesh::MeshContinuum& grid);
};

