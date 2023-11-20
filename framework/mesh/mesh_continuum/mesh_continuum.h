#pragma once

#include <memory>
#include <array>

#include "framework/mesh/mesh.h"
#include "framework/mesh/mesh_continuum/mesh_continuum_local_cell_handler.h"
#include "framework/mesh/mesh_continuum/mesh_continuum_global_cell_handler.h"
#include "framework/mesh/mesh_continuum/mesh_continuum_vertex_handler.h"

#include "framework/mpi/mpi.h"

namespace opensn
{
template <typename T>
class NDArray;
class MPICommunicatorSet;
class GridFaceHistogram;
class MeshGenerator;

/**
 * Stores the relevant information for completely defining a computationaldomain.
 */
class MeshContinuum
{
private:
  typedef std::shared_ptr<MPICommunicatorSet> MPILocalCommSetPtr;

private:
  std::vector<std::unique_ptr<Cell>> local_cells_; ///< Actual local cells
  std::vector<std::unique_ptr<Cell>> ghost_cells_; ///< Locally stored ghosts

  std::map<uint64_t, uint64_t> global_cell_id_to_local_id_map_;
  std::map<uint64_t, uint64_t> global_cell_id_to_nonlocal_id_map_;

  uint64_t global_vertex_count_ = 0;

public:
  VertexHandler vertices;
  LocalCellHandler local_cells;
  GlobalCellHandler cells;

private:
  MeshAttributes attributes = NONE;

  struct
  {
    size_t Nx = 0;
    size_t Ny = 0;
    size_t Nz = 0;
  } ortho_attributes;

  std::map<uint64_t, std::string> boundary_id_map_;

public:
  MeshContinuum()
    : local_cells(local_cells_),
      cells(local_cells_,
            ghost_cells_,
            global_cell_id_to_local_id_map_,
            global_cell_id_to_nonlocal_id_map_)
  {
  }

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

  /**Method to be called if cells and nodes have been transferred
   * to another grid.*/
  void ClearCellReferences()
  {
    local_cells_.clear();
    ghost_cells_.clear();
    global_cell_id_to_local_id_map_.clear();
    global_cell_id_to_nonlocal_id_map_.clear();
    vertices.Clear();
  }

  /**Export cells to python.
   *
   * \todo Export Cells to OBJ needs polygon support.
   */
  void ExportCellsToObj(const char* fileName, bool per_material = false, int options = 0) const;

  /**
   * Exports just the mesh to VTK format.
   */
  void ExportCellsToVTK(const std::string& file_base_name) const;

  /**
   * Exports just the portion of the mesh to ExodusII format.
   */
  void ExportCellsToExodus(const std::string& file_base_name,
                           bool suppress_node_sets = false,
                           bool suppress_side_sets = false) const;

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

  /**
   * Creates a mapping of the current face local-ids to the adjacent face's local ids.
   */
  void FindAssociatedVertices(const CellFace& cur_face, std::vector<short>& dof_mapping) const;

  /**
   * Creates a mapping of the current face local-ids to the adjacent cell's local ids.
   */
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

  /**
   * Computes the centroid from nodes specified by the given list.
   */
  Vector3 ComputeCentroidFromListOfNodes(const std::vector<uint64_t>& list) const;

  /**
   * Gets the communicator-set for interprocess communication,
   * associated with this mesh. If not created yet, it will create it.
   */
  MPILocalCommSetPtr MakeMPILocalCommunicatorSet() const;

  /**
   * Returns the total number of global cells.
   */
  size_t GetGlobalNumberOfCells() const;

  /**
   * Builds and returns a vector of unique boundary id's present in
   * the mesh.
   */
  std::vector<uint64_t> GetDomainUniqueBoundaryIDs() const;

  /**
   * Counts the number of cells within a logical volume across all partitions.
   */
  size_t CountCellsInLogicalVolume(const LogicalVolume& log_vol) const;

  /**
   * Checks whether a point is within a cell.
   */
  bool CheckPointInsideCell(const Cell& cell, const Vector3& point) const;

  MeshAttributes Attributes() const { return attributes; }

  /**
   * Gets and orthogonal mesh interface object.
   */
  std::array<size_t, 3> GetIJKInfo() const;

  /**
   * Provides a mapping from cell ijk indices to global ids.
   */
  NDArray<uint64_t> MakeIJKToGlobalIDMapping() const;

  /**
   * Determines the bounding box size of each cell and returns it as
   * a list of 3-component vectors, one Vec3 for each cell.
   */
  std::vector<Vector3> MakeCellOrthoSizes() const;

  std::pair<Vector3, Vector3> GetLocalBoundingBox() const;

private:
  friend class VolumeMesher;
  friend class MeshGenerator;
  void SetAttributes(MeshAttributes new_attribs, std::array<size_t, 3> ortho_Nis = {0, 0, 0})
  {
    attributes = attributes | new_attribs;
    ortho_attributes = {ortho_Nis[0], ortho_Nis[1], ortho_Nis[2]};
  }
};

} // namespace opensn
