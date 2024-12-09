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

/// Encapsulates all the necessary information required to fully define a computational domain.
class MeshContinuum : public std::enable_shared_from_this<MeshContinuum>
{
public:
  MeshContinuum();

  unsigned int GetDimension() const { return dim_; }
  void SetDimension(unsigned int dim) { dim_ = dim; }

  void SetGlobalVertexCount(const uint64_t count) { global_vertex_count_ = count; }
  uint64_t GetGlobalVertexCount() const { return global_vertex_count_; }

  std::map<uint64_t, std::string>& GetBoundaryIDMap() { return boundary_id_map_; }

  const std::map<uint64_t, std::string>& GetBoundaryIDMap() const { return boundary_id_map_; }

  /**
   * Makes a bndry id given a name. If the bndry name already exists,
   * the associated bndry id will be returned. Other the id will be set
   * to one more than the maximum boundary id.
   */
  uint64_t MakeBoundaryID(const std::string& boundary_name) const;

  static std::shared_ptr<MeshContinuum> New() { return std::make_shared<MeshContinuum>(); }

  /// Method to be called if cells and nodes have been transferred to another grid.
  void ClearCellReferences()
  {
    local_cells_.clear();
    ghost_cells_.clear();
    global_cell_id_to_local_id_map_.clear();
    global_cell_id_to_nonlocal_id_map_.clear();
    vertices.Clear();
  }

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

  /**
   * Check whether a cell is local by attempting to find the key in
   * the native index map.
   */
  bool IsCellLocal(uint64_t cell_global_index) const;

  /**
   * Check whether a cell is a boundary by checking if the key is found in the native or foreign
   * cell maps.
   */
  static int GetCellDimension(const Cell& cell);

  /// Creates a mapping of the current face local-ids to the adjacent face's local ids.
  void FindAssociatedVertices(const CellFace& cur_face, std::vector<short>& dof_mapping) const;

  /// Creates a mapping of the current face local-ids to the adjacent cell's local ids.
  void FindAssociatedCellVertices(const CellFace& cur_face, std::vector<short>& dof_mapping) const;

  /**
   * Given the current cell, cell A, and its adjacent cell, cell B, with
   * cell B adjacent to A at the `f`-th face of cell A. Will determine the
   * `af`-th index of the face on cell B that interface with the `f`-th face
   * of cell A.
   */
  static size_t MapCellFace(const Cell& cur_cell, const Cell& adj_cell, unsigned int f);

  /**
   * Given a global-id of a cell, will return the local-id if the cell is local, otherwise will
   * throw out_of_range.
   */
  size_t MapCellGlobalID2LocalID(uint64_t global_id) const;

  /// Computes the centroid from nodes specified by the given list.
  Vector3 ComputeCentroidFromListOfNodes(const std::vector<uint64_t>& list) const;

  /**
   * Gets the communicator-set for interprocess communication,
   * associated with this mesh. If not created yet, it will create it.
   */
  std::shared_ptr<MPICommunicatorSet> MakeMPILocalCommunicatorSet() const;

  /// Returns the total number of global cells.
  size_t GetGlobalNumberOfCells() const;

  /// Builds and returns a vector of unique boundary id's present in the mesh.
  std::vector<uint64_t> GetDomainUniqueBoundaryIDs() const;

  /**
   * Counts the number of cells within a logical volume across all partitions.
   */
  size_t CountCellsInLogicalVolume(const LogicalVolume& log_vol) const;

  /// Checks whether a point is within a cell.
  bool CheckPointInsideCell(const Cell& cell, const Vector3& point) const;

  MeshType GetType() const { return mesh_type_; }

  void SetType(MeshType type) { mesh_type_ = type; }

  bool Extruded() const { return extruded_; }

  void SetExtruded(bool extruded) { extruded_ = extruded; }

  /// Gets and orthogonal mesh interface object.
  std::array<size_t, 3> GetIJKInfo() const;

  /// Provides a mapping from cell ijk indices to global ids.
  NDArray<uint64_t, 3> MakeIJKToGlobalIDMapping() const;

  /**
   * Determines the bounding box size of each cell and returns it as
   * a list of 3-component vectors, one Vector3 for each cell.
   */
  std::vector<Vector3> MakeCellOrthoSizes() const;

  std::pair<Vector3, Vector3> GetLocalBoundingBox() const;

  /// Sets material id's for all cells to the specified material id.
  void SetUniformMaterialID(int mat_id);

  /// Sets material id's using a logical volume.
  void SetMaterialIDFromLogical(const LogicalVolume& log_vol, int mat_id, bool sense);

  /// Sets boundary id's using a logical volume.
  void SetBoundaryIDFromLogical(const LogicalVolume& log_vol,
                                const std::string& boundary_name,
                                bool sense = true);

  void SetOrthoAttributes(const OrthoMeshAttributes& attrs) { ortho_attributes_ = attrs; }

  /// Compute volume per material id's
  void ComputeVolumePerMaterialID() const;

  void SetupOrthogonalBoundaries();

private:
  /// Spatial dimension
  unsigned int dim_;
  MeshType mesh_type_;
  bool extruded_;
  OrthoMeshAttributes ortho_attributes_;
  std::map<uint64_t, std::string> boundary_id_map_;

  uint64_t global_vertex_count_;

  std::vector<std::shared_ptr<Cell>> local_cells_; ///< Actual local cells
  std::vector<std::shared_ptr<Cell>> ghost_cells_; ///< Locally stored ghosts

  std::map<uint64_t, uint64_t> global_cell_id_to_local_id_map_;
  std::map<uint64_t, uint64_t> global_cell_id_to_nonlocal_id_map_;

public:
  VertexHandler vertices;
  LocalCellHandler local_cells;
  GlobalCellHandler cells;
};

} // namespace opensn
