// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include "framework/mesh/mesh/cell.h"
#include "framework/math/geometry.h"
#include <map>
#include <array>
#include <set>

namespace opensn
{

enum MeshType : int
{
  ORTHOGONAL,
  UNSTRUCTURED
};

struct OrthoMeshAttributes
{
  size_t Nx = 0;
  size_t Ny = 0;
  size_t Nz = 0;
};

/// This object is intended for unpartitioned meshes that still require partitioning.
class UnpartitionedMesh
{
public:
  struct LightWeightFace
  {
    std::vector<uint64_t> vertex_ids;
    bool has_neighbor = false;
    uint64_t neighbor = std::numeric_limits<uint64_t>::max();

    LightWeightFace() = default;
    explicit LightWeightFace(std::vector<uint64_t> vertex_ids) : vertex_ids(std::move(vertex_ids))
    {
    }
  };
  struct LightWeightCell
  {
    const CellType type;
    const CellType sub_type;
    Vector3 centroid;
    unsigned int block_id = std::numeric_limits<unsigned int>::max();
    std::vector<uint64_t> vertex_ids;
    std::vector<LightWeightFace> faces;

    explicit LightWeightCell(CellType type, CellType sub_type) : type(type), sub_type(sub_type) {}
  };

  struct Options
  {
    std::string file_name;
    std::string block_id_fieldname = "BlockID";
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

  unsigned int GetDimension() const { return dim_; }
  void SetDimension(unsigned int dim) { dim_ = dim; }

  CoordinateSystemType GetCoordinateSystem() const { return coord_sys_; }
  void SetCoordinateSystem(const CoordinateSystemType coord_sys) { coord_sys_ = coord_sys; }

  void SetType(MeshType type) { mesh_type_ = type; }
  const MeshType& GetType() const { return mesh_type_; }

  void SetExtruded(bool extruded) { extruded_ = extruded; }
  bool IsExtruded() const { return extruded_; }

  void SetOrthoAttributes(size_t nx, size_t ny, size_t nz);
  const OrthoMeshAttributes& GetOrthoAttributes() const { return ortho_attrs_; }

  const BoundBox& GetBoundingBox() const { return bound_box_; }
  void ComputeBoundingBox();

  const std::map<uint64_t, std::string>& GetBoundaryIDMap() const { return boundary_id_map_; }
  std::map<uint64_t, std::string>& GetBoundaryIDMap() { return boundary_id_map_; }

  const std::map<std::string, uint64_t>& GetBoundaryNameMap() const { return boundary_name_map_; }
  std::map<std::string, uint64_t>& GetBoundaryNameMap() { return boundary_name_map_; }

  const std::vector<std::set<uint64_t>>& GetVertextCellSubscriptions() const
  {
    return vertex_cell_subscriptions_;
  }

  void AddCell(const std::shared_ptr<LightWeightCell>& cell) { cells_.push_back(cell); }
  size_t GetNumberOfCells() const { return cells_.size(); }

  std::vector<std::shared_ptr<LightWeightCell>>& GetCells() { return cells_; }
  const std::vector<std::shared_ptr<LightWeightCell>>& GetCells() const { return cells_; }

  std::vector<std::shared_ptr<LightWeightCell>>& GetBoundaryCells() { return boundary_cells_; }
  const std::vector<std::shared_ptr<LightWeightCell>>& GetBoundaryCells() const
  {
    return boundary_cells_;
  }

  const std::vector<Vector3>& GetVertices() const { return vertices_; }
  std::vector<Vector3>& GetVertices() { return vertices_; }

  /**
   * Establishes neighbor connectivity for the light-weight mesh.
   *
   * Two faces are connected when they have the same vertices. Throws std::logic_error when the
   * mesh is folded, that is, when the two cells of an interior face lie on the same side of it and
   * overlap, or when a face belongs to more than two cells.
   */
  void BuildMeshConnectivity();

  /// Compute centroids for all cells.
  void ComputeCentroids();

  /**
   * Checks the quality of the cells.
   *
   * Throws std::logic_error when a polygon or polyhedron is inverted, with its faces ordered
   * inward, or has zero area or volume, or when an edge is shorter than 1e-12. Warns when the
   * centroid of a cell lies outside the plane of a triangle of one of its faces, as it can in a
   * non-convex cell.
   */
  void CheckQuality();

  /// Makes or gets a boundary that uniquely identifies the given name.
  uint64_t MakeBoundaryID(const std::string& boundary_name);

  void AddBoundary(uint64_t id, const std::string& name);

protected:
  /// Spatial mesh dimension
  unsigned int dim_;
  CoordinateSystemType coord_sys_;
  MeshType mesh_type_;
  bool extruded_;
  OrthoMeshAttributes ortho_attrs_;
  Options mesh_options_;
  BoundBox bound_box_;
  std::map<uint64_t, std::string> boundary_id_map_;
  std::map<std::string, uint64_t> boundary_name_map_;

  std::vector<Vector3> vertices_;
  std::vector<std::shared_ptr<LightWeightCell>> cells_;
  std::vector<std::shared_ptr<LightWeightCell>> boundary_cells_;
  std::vector<std::set<uint64_t>> vertex_cell_subscriptions_;
};

} // namespace opensn
