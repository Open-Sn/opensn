#pragma once

#include <memory>
#include <array>

#include "../chi_mesh.h"
#include "chi_meshcontinuum_localcellhandler.h"
#include "chi_meshcontinuum_globalcellhandler.h"
#include "chi_meshcontinuum_vertexhandler.h"

#include "chi_mpi.h"

namespace chi_data_types
{
template <typename T>
class NDArray;
} // namespace chi_data_types

namespace chi
{
class ChiMPICommunicatorSet;
}

namespace chi_mesh
{
class GridFaceHistogram;
class MeshGenerator;
} // namespace chi_mesh

// ######################################################### Class Definition
/**Stores the relevant information for completely defining a computational
 * domain. */
class chi_mesh::MeshContinuum
{
private:
  typedef std::shared_ptr<chi::ChiMPICommunicatorSet> MPILocalCommSetPtr;

private:
  std::vector<std::unique_ptr<chi_mesh::Cell>> local_cells_; ///< Actual local cells
  std::vector<std::unique_ptr<chi_mesh::Cell>> ghost_cells_; ///< Locally stored ghosts

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

  void ExportCellsToObj(const char* fileName, bool per_material = false, int options = 0) const;
  void ExportCellsToVTK(const std::string& file_base_name) const;
  void ExportCellsToExodus(const std::string& file_base_name,
                           bool suppress_node_sets = false,
                           bool suppress_side_sets = false) const;

  std::shared_ptr<GridFaceHistogram> MakeGridFaceHistogram(double master_tolerance = 100.0,
                                                           double slave_tolerance = 1.1) const;

  bool IsCellLocal(uint64_t cell_global_index) const;
  static int GetCellDimension(const chi_mesh::Cell& cell);

  /**Creates a mapping of the current face local-ids to the
   * adjacent face's local ids.*/
  void FindAssociatedVertices(const chi_mesh::CellFace& cur_face,
                              std::vector<short>& dof_mapping) const;
  /**Creates a mapping of the current face local-ids to the
   * adjacent cell's local ids.*/
  void FindAssociatedCellVertices(const chi_mesh::CellFace& cur_face,
                                  std::vector<short>& dof_mapping) const;
  static size_t
  MapCellFace(const chi_mesh::Cell& cur_cell, const chi_mesh::Cell& adj_cell, unsigned int f);

  /**Given a global-id of a cell, will return the local-id if the
   * cell is local, otherwise will throw out_of_range.*/
  size_t MapCellGlobalID2LocalID(uint64_t global_id) const;

  chi_mesh::Vector3 ComputeCentroidFromListOfNodes(const std::vector<uint64_t>& list) const;

  MPILocalCommSetPtr MakeMPILocalCommunicatorSet() const;

  size_t GetGlobalNumberOfCells() const;

  std::vector<uint64_t> GetDomainUniqueBoundaryIDs() const;

  size_t CountCellsInLogicalVolume(const chi_mesh::LogicalVolume& log_vol) const;
  bool CheckPointInsideCell(const chi_mesh::Cell& cell, const chi_mesh::Vector3& point) const;

  MeshAttributes Attributes() const { return attributes; }

  std::array<size_t, 3> GetIJKInfo() const;
  chi_data_types::NDArray<uint64_t> MakeIJKToGlobalIDMapping() const;
  std::vector<chi_mesh::Vector3> MakeCellOrthoSizes() const;

  std::pair<chi_mesh::Vector3, chi_mesh::Vector3> GetLocalBoundingBox() const;

private:
  friend class chi_mesh::VolumeMesher;
  friend class chi_mesh::MeshGenerator;
  void SetAttributes(MeshAttributes new_attribs, std::array<size_t, 3> ortho_Nis = {0, 0, 0})
  {
    attributes = attributes | new_attribs;
    ortho_attributes = {ortho_Nis[0], ortho_Nis[1], ortho_Nis[2]};
  }
};

