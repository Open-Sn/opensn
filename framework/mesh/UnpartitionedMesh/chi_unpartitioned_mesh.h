#ifndef CHI_MESH_UNPARTITIONED_MESH_H
#define CHI_MESH_UNPARTITIONED_MESH_H

#include "mesh/chi_mesh.h"
#include "mesh/Cell/cell.h"

class vtkCell;
class vtkUnstructuredGrid;
template <class T>
class vtkSmartPointer;

#include <map>
#include <array>

// ###################################################################
/**This object is intented for unpartitioned meshes that still require
 * partitioning.*/
class chi_mesh::UnpartitionedMesh
{
public:
  struct LightWeightFace
  {
    std::vector<uint64_t> vertex_ids;
    bool has_neighbor = false;
    uint64_t neighbor = 0;

    LightWeightFace() = default;
    explicit LightWeightFace(std::vector<uint64_t> in_vertex_ids)
      : vertex_ids(std::move(in_vertex_ids))
    {
    }
  };
  struct LightWeightCell
  {
    const chi_mesh::CellType type;
    const chi_mesh::CellType sub_type;
    chi_mesh::Vertex centroid;
    int material_id = -1;
    std::vector<uint64_t> vertex_ids;
    std::vector<LightWeightFace> faces;

    explicit LightWeightCell(chi_mesh::CellType in_type,
                             chi_mesh::CellType in_sub_type)
      : type(in_type), sub_type(in_sub_type)
    {
    }
  };

  struct Options
  {
    std::string file_name;
    std::string material_id_fieldname = "BlockID";
    std::string boundary_id_fieldname;
    double scale = 1.0;
    size_t ortho_Nx = 0;
    size_t ortho_Ny = 0;
    size_t ortho_Nz = 0;

    std::map<uint64_t, std::string> boundary_id_map;
  };

  struct BoundBox
  {
    double xmin = 0.0, xmax = 0.0, ymin = 0.0, ymax = 0.0, zmin = 0.0,
           zmax = 0.0;
  };

protected:
  std::vector<chi_mesh::Vertex> vertices_;
  std::vector<LightWeightCell*> raw_cells_;
  std::vector<LightWeightCell*> raw_boundary_cells_;
  std::vector<std::set<uint64_t>> vertex_cell_subscriptions_;

  MeshAttributes attributes_ = NONE;
  Options mesh_options_;
  std::shared_ptr<BoundBox> bound_box_ = nullptr;

protected:
  static LightWeightCell* CreateCellFromVTKPolyhedron(vtkCell* vtk_cell);
  static LightWeightCell* CreateCellFromVTKPolygon(vtkCell* vtk_cell);
  static LightWeightCell* CreateCellFromVTKLine(vtkCell* vtk_cell);
  static LightWeightCell* CreateCellFromVTKVertex(vtkCell* vtk_cell);

  typedef vtkSmartPointer<vtkUnstructuredGrid> vtkUGridPtr;
  typedef std::pair<vtkUGridPtr, std::string> vtkUGridPtrAndName;
  void CopyUGridCellsAndPoints(vtkUnstructuredGrid& ugrid,
                               double scale,
                               int dimension_to_copy);

  void SetMaterialIDsFromList(const std::vector<int>& material_ids);

  void
  SetBoundaryIDsFromBlocks(std::vector<vtkUGridPtrAndName>& bndry_grid_blocks);

public:
  const BoundBox& GetBoundBox() const { return *bound_box_; }

  Options& GetMeshOptions() { return mesh_options_; }
  const Options& GetMeshOptions() const { return mesh_options_; }

  MeshAttributes& GetMeshAttributes() { return attributes_; }
  const MeshAttributes& GetMeshAttributes() const { return attributes_; }

  const std::vector<std::set<uint64_t>>& GetVertextCellSubscriptions() const
  {
    return vertex_cell_subscriptions_;
  }

  void AddCell(LightWeightCell*& cell) { raw_cells_.push_back(cell); }
  size_t GetNumberOfCells() const { return raw_cells_.size(); }
  std::vector<LightWeightCell*>& GetRawCells() { return raw_cells_; }
  const std::vector<LightWeightCell*>& GetRawCells() const
  {
    return raw_cells_;
  }

  const std::vector<chi_mesh::Vertex>& GetVertices() const { return vertices_; }
  std::vector<chi_mesh::Vertex>& GetVertices() { return vertices_; }

  void BuildMeshConnectivity();
  void ComputeCentroidsAndCheckQuality();
  /**Makes or gets a boundary that uniquely identifies the given name.*/
  uint64_t MakeBoundaryID(const std::string& boundary_name);
  void SetAttributes(MeshAttributes new_attribs,
                     std::array<size_t, 3> ortho_Nis = {0, 0, 0});

  void ReadFromVTU(const Options& options);
  void ReadFromPVTU(const Options& options);
  void ReadFromEnsightGold(const Options& options);
  void ReadFromWavefrontOBJ(const Options& options);

  void ReadFromMsh(const Options& options);

  void ReadFromExodus(const Options& options);

  void PushProxyCell(const std::string& type_str,
                     const std::string& sub_type_str,
                     int cell_num_faces,
                     int cell_material_id,
                     const std::vector<std::vector<uint64_t>>& proxy_faces);

  ~UnpartitionedMesh();

  void CleanUp()
  {
    for (auto& cell : raw_cells_)
      delete cell;
    for (auto& cell : raw_boundary_cells_)
      delete cell;
    vertices_.clear();
    vertices_.shrink_to_fit();
    raw_cells_.clear();
    raw_cells_.shrink_to_fit();
    raw_boundary_cells_.clear();
    raw_boundary_cells_.shrink_to_fit();
    vertex_cell_subscriptions_.clear();
    vertex_cell_subscriptions_.shrink_to_fit();
  }
};

#endif // CHI_MESH_UNPARTITIONED_MESH_H
