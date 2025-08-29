// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include "framework/data_types/ndarray.h"
#include "framework/mesh/mesh.h"
#include "framework/mesh/mesh_continuum/mesh_continuum_local_cell_handler.h"
#include "framework/mesh/mesh_continuum/mesh_continuum_global_cell_handler.h"
#include "framework/mesh/mesh_continuum/mesh_continuum_vertex_handler.h"
#include <memory>
#include <array>

namespace opensn
{
class MPICommunicatorSet;
class GridFaceHistogram;
class MeshGenerator;
class LogicalVolume;

/// Encapsulates all the necessary information required to fully define a computational domain.
class MeshContinuum
{
public:
  MeshContinuum();

  unsigned int GetDimension() const { return dim_; }
  void SetDimension(const unsigned int dim) { dim_ = dim; }

  MeshType GetType() const { return mesh_type_; }
  void SetType(const MeshType type) { mesh_type_ = type; }

  CoordinateSystemType GetCoordinateSystem() const { return coord_sys_; }
  void SetCoordinateSystem(const CoordinateSystemType coord_sys) { coord_sys_ = coord_sys; }

  bool Extruded() const { return extruded_; }
  void SetExtruded(const bool extruded) { extruded_ = extruded; }

  /// Gets and orthogonal mesh interface object.
  std::array<size_t, 3> GetIJKInfo() const;
  void SetOrthoAttributes(const OrthoMeshAttributes& attrs) { ortho_attributes_ = attrs; }

  void SetGlobalVertexCount(const uint64_t count) { global_vertex_count_ = count; }
  uint64_t GetGlobalVertexCount() const { return global_vertex_count_; }
  size_t GetGlobalNumberOfCells() const;

  std::map<uint64_t, std::string>& GetBoundaryIDMap() { return boundary_id_map_; }
  const std::map<uint64_t, std::string>& GetBoundaryIDMap() const { return boundary_id_map_; }
  /// Returns the unique boundary ids present in the problem.
  std::vector<uint64_t> GetUniqueBoundaryIDs() const;

  /// Compute the geometric data for the cells and faces in the mesh.
  void ComputeGeometricInfo();

  /// Method to be called if cells and nodes have been transferred to another grid.
  void ClearCellReferences();

  /**
   * Makes a boundary id given a name. If the boundary name already exists, the associated
   * boundary id will be returned. Other the id will be set to one more than the maximum boundary
   * id.
   */
  uint64_t MakeBoundaryID(const std::string& boundary_name) const;

  /// Defines the standard x/y/z min/max boundaries.
  void SetOrthogonalBoundaries();

  /**
   * Populates a face histogram.
   *
   * \param master_tolerance Multiple histograms will only be attempted
   * if the ratio of the maximum dofs-per-face to the average dofs-per-face
   * is greater than this value. Default 1.2.
   *
   * \param slave_tolerance While traversing a sorted list of dofs-per-face,
   * a new bin will only be generated when the ratio of the listed dofs-per-face
   * to a running bin average exceeds this value. Defualt 1.1.
   *
   * The function populates face_categories which is a structure containing
   * pairs. Pair.first is the max dofs-per-face for the category and Pair.second
   * is the number of faces in this category.
   */
  std::shared_ptr<GridFaceHistogram> MakeGridFaceHistogram(double master_tolerance = 100.0,
                                                           double slave_tolerance = 1.1) const;

  /// Returns whether the cell with the given global id is locally owned.
  bool IsCellLocal(uint64_t global_id) const;
  /**
   * Given a global id of a cell, returns the local id if the cell is local, otherwise throws
   * an out_of_range error.
   */
  size_t MapCellGlobalID2LocalID(uint64_t global_id) const;

  /// Creates a mapping of the current face local ids to the adjacent face's local ids.
  void FindAssociatedVertices(const CellFace& cur_face, std::vector<short>& dof_mapping) const;
  /// Creates a mapping of the current face local ids to the adjacent cell's local ids.
  void FindAssociatedCellVertices(const CellFace& cur_face, std::vector<short>& dof_mapping) const;

  /// Counts the number of cells within a logical volume across all partitions.
  size_t CountCellsInLogicalVolume(const LogicalVolume& log_vol) const;

  /// Checks whether a point is within a cell.
  bool CheckPointInsideCell(const Cell& cell, const Vector3& point) const;
  /// Checks whether a point is within a cell face.
  bool CheckPointInsideCellFace(const Cell& cell, size_t face_i, const Vector3& point) const;

  /// Provides a mapping from cell ijk indices to global ids.
  NDArray<uint64_t, 3> MakeIJKToGlobalIDMapping() const;

  /**
   * Determines the bounding box size of each cell and returns it as a list of 3-component vectors,
   * one Vector3 for each cell.
   */
  std::vector<Vector3> MakeCellOrthoSizes() const;

  /// Returns the bounding box corners for the locally owned cells.
  std::pair<Vector3, Vector3> GetLocalBoundingBox() const;

  /// Sets block ids for all cells to the specified block id.
  void SetUniformBlockID(int blk_id);

  /// Sets block IDs using a logical volume.
  void SetBlockIDFromLogicalVolume(const LogicalVolume& log_vol, int blk_id, bool sense);

  /// Sets boundary ids using a logical volume.
  void SetBoundaryIDFromLogicalVolume(const LogicalVolume& log_vol,
                                      const std::string& boundary_name,
                                      bool sense = true);

  /// Computes the centroid from nodes specified by the given list.
  Vector3 ComputeCentroidFromListOfNodes(const std::vector<uint64_t>& list) const;

  /**
   * Get the face vertices of a tetrahedron contained within the given face and
   * side of a polyhedron.
   */
  std::array<std::array<Vector3, 3>, 4>
  GetTetrahedralFaceVertices(const Cell& cell, const CellFace& face, size_t side) const;

  /**
   * Gets the communicator-set for interprocess communication, associated with this mesh.
   * If not created yet, it will create it.
   */
  std::shared_ptr<MPICommunicatorSet> MakeMPILocalCommunicatorSet() const;

  VertexHandler vertices;
  LocalCellHandler local_cells;
  GlobalCellHandler cells;

  /// Compute volume per block IDs
  std::map<int, double> ComputeVolumePerBlockID() const;

private:
  /// Spatial dimension
  unsigned int dim_;
  MeshType mesh_type_;
  CoordinateSystemType coord_sys_;
  bool extruded_;
  OrthoMeshAttributes ortho_attributes_;
  std::map<uint64_t, std::string> boundary_id_map_;

  uint64_t global_vertex_count_;

  /// Locally owned cells
  std::vector<std::shared_ptr<Cell>> local_cells_;
  /// Locally stored ghost cells
  std::vector<std::shared_ptr<Cell>> ghost_cells_;

  std::map<uint64_t, uint64_t> global_cell_id_to_local_id_map_;
  std::map<uint64_t, uint64_t> global_cell_id_to_nonlocal_id_map_;

public:
  /// Returns a new instance of the spatial discretization.
  static std::shared_ptr<MeshContinuum> New() { return std::make_shared<MeshContinuum>(); }

  /// Returns the spatial dimensionality of the cell.
  static int GetCellDimension(const Cell& cell);

  /**
   * Given the current cell, cell A, and its adjacent cell, cell B, with cell B adjacent to
   * A at the `f`-th face of cell A. Will determine the `af`-th index of the face on cell B
   * that interface with the `f`-th face of cell A.
   */
  static size_t MapCellFace(const Cell& cur_cell, const Cell& adj_cell, unsigned int f);
};

} // namespace opensn
