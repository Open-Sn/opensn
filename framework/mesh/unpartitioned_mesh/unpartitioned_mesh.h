// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include "framework/mesh/cell/cell.h"
#include "framework/mesh/mesh.h"
#include <map>
#include <array>
#include <set>

namespace opensn
{

/**
 * This object is intended for unpartitioned meshes that still require partitioning.
 */
class UnpartitionedMesh
{
public:
  struct LightWeightFace
  {
    std::vector<uint64_t> vertex_ids;
    bool has_neighbor = false;
    uint64_t neighbor = 0;

    LightWeightFace() = default;
    explicit LightWeightFace(std::vector<uint64_t> vertex_ids) : vertex_ids(std::move(vertex_ids))
    {
    }
  };
  struct LightWeightCell
  {
    const CellType type;
    const CellType sub_type;
    Vertex centroid;
    int material_id = -1;
    std::vector<uint64_t> vertex_ids;
    std::vector<LightWeightFace> faces;

    explicit LightWeightCell(CellType type, CellType sub_type) : type(type), sub_type(sub_type) {}
  };

  struct Options
  {
    std::string file_name;
    std::string material_id_fieldname = "BlockID";
    std::string boundary_id_fieldname;
    double scale = 1.0;
  };

  struct BoundBox
  {
    double xmin = 0.0, xmax = 0.0, ymin = 0.0, ymax = 0.0, zmin = 0.0, zmax = 0.0;
  };

public:
  UnpartitionedMesh();
  ~UnpartitionedMesh();

  unsigned int Dimension() const { return dim_; }
  void SetDimension(unsigned int dim) { dim_ = dim; }

  const BoundBox& BoundingBox() const { return bound_box_; }
  void ComputeBoundingBox();

  void SetMeshOptions(const Options& opts) { mesh_options_ = opts; }
  const Options& MeshOptions() const { return mesh_options_; }

  void SetAttributes(MeshAttributes attrs) { attributes_ = attrs; }
  const MeshAttributes& Attributes() const { return attributes_; }

  const std::vector<std::set<uint64_t>>& GetVertextCellSubscriptions() const
  {
    return vertex_cell_subscriptions_;
  }

  void AddCell(LightWeightCell*& cell) { raw_cells_.push_back(cell); }
  size_t GetNumberOfCells() const { return raw_cells_.size(); }

  std::vector<LightWeightCell*>& RawCells() { return raw_cells_; }
  const std::vector<LightWeightCell*>& RawCells() const { return raw_cells_; }

  std::vector<LightWeightCell*>& RawBoundaryCells() { return raw_boundary_cells_; }
  const std::vector<LightWeightCell*>& RawBoundaryCells() const { return raw_boundary_cells_; }

  const std::vector<Vertex>& Vertices() const { return vertices_; }
  std::vector<Vertex>& Vertices() { return vertices_; }

  /**
   * Establishes neighbor connectivity for the light-weight mesh.
   */
  void BuildMeshConnectivity();

  /**
   * Compute centroids for all cells.
   */
  void ComputeCentroids();

  /**
   * Check element quality
   */
  void CheckQuality();

  /**
   * Makes or gets a boundary that uniquely identifies the given name.
   */
  uint64_t MakeBoundaryID(const std::string& boundary_name);

  void AddBoundary(uint64_t id, const std::string& name);

  const std::map<uint64_t, std::string>& BoundaryIDMap() const { return boundary_id_map_; }
  std::map<uint64_t, std::string>& BoundaryIDMap() { return boundary_id_map_; }

  void SetOrthoAttributes(size_t nx, size_t ny, size_t nz);
  const OrthoMeshAttributes& OrthoAttributes() const { return ortho_attrs_; }

  /**Makes a cell from proxy information and pushes the cell to the mesh.*/
  void PushProxyCell(const std::string& type_str,
                     const std::string& sub_type_str,
                     int cell_num_faces,
                     int cell_material_id,
                     const std::vector<std::vector<uint64_t>>& proxy_faces);

  void CleanUp();

protected:
  /// Spatial mesh dimension
  unsigned int dim_;
  MeshAttributes attributes_;
  OrthoMeshAttributes ortho_attrs_;
  Options mesh_options_;
  BoundBox bound_box_;
  std::map<uint64_t, std::string> boundary_id_map_;

  std::vector<Vertex> vertices_;
  std::vector<LightWeightCell*> raw_cells_;
  std::vector<LightWeightCell*> raw_boundary_cells_;
  std::vector<std::set<uint64_t>> vertex_cell_subscriptions_;
};

} // namespace opensn
