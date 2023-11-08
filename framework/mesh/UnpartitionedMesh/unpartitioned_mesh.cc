#include "framework/mesh/UnpartitionedMesh/unpartitioned_mesh.h"
#include "framework/mesh/Cell/cell.h"
#include "framework/mesh/MeshContinuum/chi_grid_vtk_utils.h"
#include "framework/chi_runtime.h"
#include "framework/logging/chi_log.h"
#include "framework/utils/chi_timer.h"
#include "framework/mpi/chi_mpi.h"
#include "framework/utils/chi_utils.h"
#include <vtkPolygon.h>
#include <vtkLine.h>
#include <vtkVertex.h>
#include <vtkUnstructuredGrid.h>
#include <vtkCellData.h>
#include <vtkPointData.h>
#include <vtkSmartPointer.h>
#include <vtkXMLUnstructuredGridReader.h>
#include <vtkInformation.h>
#include <vtkXMLPUnstructuredGridReader.h>
#include <vtkEnSightGoldBinaryReader.h>
#include <vtkMultiBlockDataSet.h>
#include <vtkExodusIIReader.h>
#include <algorithm>
#include <fstream>
#include <map>

#define ErrorReadingFile(fname)                                                                    \
  std::runtime_error("Failed to open file: " + options.file_name + " in call to " + #fname + ".")

namespace chi_mesh
{

UnpartitionedMesh::~UnpartitionedMesh()
{
  for (auto& cell : raw_cells_)
    delete cell;
  for (auto& cell : raw_boundary_cells_)
    delete cell;
  Chi::log.Log0Verbose1() << "Unpartitioned Mesh destructor called";
}

void
UnpartitionedMesh::ComputeCentroidsAndCheckQuality()
{
  const Vector3 khat(0.0, 0.0, 1.0);

  Chi::log.Log0Verbose1() << "Computing cell-centroids.";
  for (auto cell : raw_cells_)
  {
    cell->centroid = Vertex(0.0, 0.0, 0.0);
    for (auto vid : cell->vertex_ids)
      cell->centroid += vertices_[vid];

    cell->centroid = cell->centroid / static_cast<double>(cell->vertex_ids.size());
  }
  Chi::log.Log0Verbose1() << "Done computing cell-centroids.";

  Chi::log.Log0Verbose1() << "Checking cell-center-to-face orientations";
  size_t num_negative_volume_elements = 0;
  for (auto cell : raw_cells_)
  {
    if (cell->type == CellType::POLYGON)
    {
      // Form triangles
      size_t num_verts = cell->vertex_ids.size();
      for (size_t v = 0; v < num_verts; ++v)
      {
        size_t vp1 = (v < (num_verts - 1)) ? v + 1 : 0;

        const auto& v0 = vertices_[cell->vertex_ids[v]];
        const auto& v1 = vertices_[cell->vertex_ids[vp1]];

        auto E01 = v1 - v0;
        auto n = E01.Cross(khat).Normalized();

        if (n.Dot(v0 - cell->centroid) < 0.0) ++num_negative_volume_elements;
      } // for v
    }
    else if (cell->type == CellType::POLYHEDRON)
    {
      for (auto& face : cell->faces)
      {
        if (face.vertex_ids.size() < 2)
          throw std::logic_error(std::string(__PRETTY_FUNCTION__) +
                                 ": cell-center-to-face check encountered face "
                                 "with less than 2 vertices on a face, making "
                                 "normal computation impossible.");

        // Compute centroid
        Vector3 face_centroid;
        for (uint64_t vid : face.vertex_ids)
          face_centroid += vertices_[vid];
        face_centroid /= static_cast<double>(face.vertex_ids.size());

        // Form tets for each face edge
        size_t num_face_verts = face.vertex_ids.size();
        for (size_t fv = 0; fv < face.vertex_ids.size(); ++fv)
        {
          size_t fvp1 = (fv < (num_face_verts - 1)) ? fv + 1 : 0;

          const auto& fv1 = vertices_[face.vertex_ids[fv]];
          const auto& fv2 = vertices_[face.vertex_ids[fvp1]];

          auto E0 = fv1 - face_centroid;
          auto E1 = fv2 - face_centroid;
          auto n = E0.Cross(E1).Normalized();

          if (n.Dot(fv1 - cell->centroid) < 0.0) ++num_negative_volume_elements;
        }
      } // for face
    }   // if polyhedron
  }     // for cell in raw_cells

  Chi::log.Log0Verbose1() << "Checking face sizes";
  size_t cell_id = 0;
  for (auto cell : raw_cells_)
  {
    if (cell->type == CellType::POLYGON)
    {
      size_t f = 0;
      for (const auto& face : cell->faces)
      {
        const auto& v0 = vertices_.at(face.vertex_ids[0]);
        const auto& v1 = vertices_.at(face.vertex_ids[1]);
        ChiLogicalErrorIf((v1 - v0).Norm() < 1.0e-12,
                          "Cell " + std::to_string(cell_id) +
                            " (centroid=" + cell->centroid.PrintStr() + ") face " +
                            std::to_string(f) + ": Face has length < 1.0e-12.");
        ++f;
      }
    } // if polygon
    if (cell->type == CellType::POLYHEDRON)
    {
      size_t f = 0;
      for (const auto& face : cell->faces)
      {
        size_t num_face_verts = face.vertex_ids.size();
        for (size_t s = 0; s < face.vertex_ids.size(); ++s)
        {
          size_t fvp1 = (s < (num_face_verts - 1)) ? s + 1 : 0;

          const auto& v0 = vertices_.at(face.vertex_ids[s]);
          const auto& v1 = vertices_.at(face.vertex_ids[fvp1]);

          ChiLogicalErrorIf((v1 - v0).Norm() < 1.0e-12,
                            "Cell " + std::to_string(cell_id) + " (centroid=" +
                              cell->centroid.PrintStr() + ") face " + std::to_string(f) + " side " +
                              std::to_string(s) + ": Side has length < 1.0e-12.");
        }

        ++f;
      }
    } // if polyhedron
    ++cell_id;
  } // for cell in raw_cells

  if (num_negative_volume_elements > 0)
    Chi::log.LogAllWarning() << "Cell quality checks detected " << num_negative_volume_elements
                             << " negative volume sub-elements (sub-triangle or sub-tetrahedron)."
                             << " This issue could result in incorrect quantities"
                             << " under some circumstances.";
  Chi::log.Log() << "Done checking cell-center-to-face orientations";
}

uint64_t
UnpartitionedMesh::MakeBoundaryID(const std::string& boundary_name)
{
  auto& boundary_id_map = mesh_options_.boundary_id_map;
  if (boundary_id_map.empty()) return 0;

  for (const auto& [id, name] : boundary_id_map)
    if (boundary_name == name) return id;

  uint64_t max_id = 0;
  for (const auto& [id, name] : boundary_id_map)
    max_id = std::max(id, max_id);

  return max_id + 1;
}

void
UnpartitionedMesh::SetAttributes(MeshAttributes new_attribs, std::array<size_t, 3> ortho_Nis)
{
  attributes_ = attributes_ | new_attribs;
  mesh_options_.ortho_Nx = ortho_Nis[0];
  mesh_options_.ortho_Ny = ortho_Nis[1];
  mesh_options_.ortho_Nz = ortho_Nis[2];
}

void
UnpartitionedMesh::BuildMeshConnectivity()
{
  const size_t num_raw_cells = raw_cells_.size();
  const size_t num_raw_vertices = vertices_.size();

  //======================================== Reset all cell neighbors
  int num_bndry_faces = 0;
  for (auto& cell : raw_cells_)
    for (auto& face : cell->faces)
      if (not face.has_neighbor) ++num_bndry_faces;

  Chi::log.Log0Verbose1() << Chi::program_timer.GetTimeString()
                          << " Number of unconnected faces "
                             "before connectivity: "
                          << num_bndry_faces;

  Chi::log.Log() << Chi::program_timer.GetTimeString() << " Establishing cell connectivity.";

  //======================================== Establish internal connectivity
  // Populate vertex subscriptions to internal cells
  vertex_cell_subscriptions_.resize(num_raw_vertices);
  {
    uint64_t cur_cell_id = 0;
    for (const auto& cell : raw_cells_)
    {
      for (auto vid : cell->vertex_ids)
        vertex_cell_subscriptions_.at(vid).insert(cur_cell_id);
      ++cur_cell_id;
    }
  }

  Chi::log.Log() << Chi::program_timer.GetTimeString() << " Vertex cell subscriptions complete.";

  // Process raw cells
  {
    uint64_t aux_counter = 0;
    uint64_t cur_cell_id = 0;
    for (auto& cell : raw_cells_)
    {
      for (auto& cur_cell_face : cell->faces)
      {
        if (cur_cell_face.has_neighbor) { continue; }
        const std::set<uint64_t> cfvids(cur_cell_face.vertex_ids.begin(),
                                        cur_cell_face.vertex_ids.end());

        std::set<size_t> cells_to_search;
        for (uint64_t vid : cfvids)
          for (uint64_t cell_id : vertex_cell_subscriptions_.at(vid))
            if (cell_id != cur_cell_id) cells_to_search.insert(cell_id);

        for (uint64_t adj_cell_id : cells_to_search)
        {
          auto adj_cell = raw_cells_.at(adj_cell_id);

          for (auto& adj_cell_face : adj_cell->faces)
          {
            if (adj_cell_face.has_neighbor) { continue; }
            const std::set<uint64_t> afvids(adj_cell_face.vertex_ids.begin(),
                                            adj_cell_face.vertex_ids.end());

            if (cfvids == afvids)
            {
              cur_cell_face.neighbor = adj_cell_id;
              adj_cell_face.neighbor = cur_cell_id;

              cur_cell_face.has_neighbor = true;
              adj_cell_face.has_neighbor = true;

              goto face_neighbor_found;
            }
          } // for adjacent cell face
        }
      face_neighbor_found:;
      } // for face

      ++cur_cell_id;
      const double fraction_complete =
        static_cast<double>(cur_cell_id) / static_cast<double>(num_raw_cells);
      if (fraction_complete >= static_cast<double>(aux_counter + 1) * 0.1)
      {
        Chi::log.Log() << Chi::program_timer.GetTimeString() << " Surpassing cell " << cur_cell_id
                       << " of " << num_raw_cells << " (" << (aux_counter + 1) * 10 << "%)";
        ++aux_counter;
      }
    } // for cell
  }

  Chi::log.Log() << Chi::program_timer.GetTimeString()
                 << " Establishing cell boundary connectivity.";

  //======================================== Establish boundary connectivity
  // Make list of internal cells on the boundary
  std::vector<LightWeightCell*> internal_cells_on_boundary;
  for (auto& cell : raw_cells_)
  {
    bool cell_on_boundary = false;
    for (auto& face : cell->faces)
      if (not face.has_neighbor)
      {
        cell_on_boundary = true;
        break;
      }

    if (cell_on_boundary) internal_cells_on_boundary.push_back(cell);
  }

  // Populate vertex subscriptions to boundary cells
  std::vector<std::set<uint64_t>> vertex_bndry_cell_subscriptions(vertices_.size());
  {
    uint64_t cur_cell_id = 0;
    for (auto& cell : raw_boundary_cells_)
    {
      for (auto vid : cell->vertex_ids)
        vertex_bndry_cell_subscriptions.at(vid).insert(cur_cell_id);
      ++cur_cell_id;
    }
  }

  // Process boundary cells
  for (auto& cell : internal_cells_on_boundary)
    for (auto& face : cell->faces)
    {
      if (face.has_neighbor) continue;
      std::set<uint64_t> cfvids(face.vertex_ids.begin(), face.vertex_ids.end());

      std::set<size_t> cells_to_search;
      for (uint64_t vid : face.vertex_ids)
        for (uint64_t cell_id : vertex_bndry_cell_subscriptions[vid])
          cells_to_search.insert(cell_id);

      for (uint64_t adj_cell_id : cells_to_search)
      {
        auto& adj_cell = raw_boundary_cells_[adj_cell_id];

        std::set<uint64_t> afvids(adj_cell->vertex_ids.begin(), adj_cell->vertex_ids.end());

        if (cfvids == afvids)
        {
          face.neighbor = adj_cell->material_id;
          break;
        }
      } // for adj_cell_id
    }   // for face

  num_bndry_faces = 0;
  for (auto cell : raw_cells_)
    for (auto& face : cell->faces)
      if (not face.has_neighbor) ++num_bndry_faces;

  Chi::log.Log0Verbose1() << Chi::program_timer.GetTimeString()
                          << " Number of boundary faces "
                             "after connectivity: "
                          << num_bndry_faces;

  Chi::log.Log() << Chi::program_timer.GetTimeString() << " Done establishing cell connectivity.";
}

UnpartitionedMesh::LightWeightCell*
UnpartitionedMesh::CreateCellFromVTKPolyhedron(vtkCell* vtk_cell)
{
  const std::string fname = "UnpartitionedMesh::CreateCellFromVTKPolyhedron";

  CellType sub_type;
  switch (vtk_cell->GetCellType())
  {
    case VTK_HEXAGONAL_PRISM:
    case VTK_PENTAGONAL_PRISM:
    case VTK_POLYHEDRON:
      sub_type = CellType::POLYHEDRON;
      break;
    case VTK_PYRAMID:
      sub_type = CellType::PYRAMID;
      break;
    case VTK_WEDGE:
      sub_type = CellType::WEDGE;
      break;
    case VTK_HEXAHEDRON:
    case VTK_VOXEL:
      sub_type = CellType::HEXAHEDRON;
      break;
    case VTK_TETRA:
      sub_type = CellType::TETRAHEDRON;
      break;
    default:
      throw std::logic_error(fname + ": Unsupported 3D cell type encountered.");
  }
  auto polyh_cell = new LightWeightCell(CellType::POLYHEDRON, sub_type);

  auto num_cpoints = vtk_cell->GetNumberOfPoints();
  auto num_cfaces = vtk_cell->GetNumberOfFaces();

  polyh_cell->vertex_ids.reserve(num_cpoints);
  auto point_ids = vtk_cell->GetPointIds();
  for (int p = 0; p < num_cpoints; ++p)
  {
    uint64_t point_id = point_ids->GetId(p);
    polyh_cell->vertex_ids.push_back(point_id);
  } // for p

  switch (sub_type)
  {
    // The cell vertex ids in VTK is the same as in ChiTech so we don't
    // need to remap the vertices. We do however need to remap the faces.
    case CellType::HEXAHEDRON:
    {
      std::vector<std::vector<uint64_t>> face_vids = {
        {1, 2, 6, 5}, {3, 0, 4, 7}, {2, 3, 7, 6}, {0, 1, 5, 4}, {4, 5, 6, 7}, {3, 2, 1, 0}};
      for (int f = 0; f < 6; ++f)
      {
        LightWeightFace face;

        face.vertex_ids.reserve(4);
        for (int p = 0; p < 4; ++p)
          face.vertex_ids.push_back(polyh_cell->vertex_ids[face_vids[f][p]]);

        polyh_cell->faces.push_back(face);
      }
      break;
    }
    // For wedges we need to remap cell vertices and faces.
    case CellType::WEDGE:
    {
      std::vector<uint64_t> remapping = {0, 2, 1, 3, 5, 4};
      std::vector<uint64_t> remapped(6, 0);
      for (int i = 0; i < 6; ++i)
        remapped[i] = polyh_cell->vertex_ids[remapping[i]];
      polyh_cell->vertex_ids = remapped;

      std::vector<std::vector<uint64_t>> face_vids = {
        {0, 1, 4, 3}, {1, 2, 5, 4}, {2, 0, 3, 5}, {3, 4, 5}, {0, 2, 1}};
      for (int f = 0; f < 6; ++f)
      {
        LightWeightFace face;

        face.vertex_ids.reserve(4);
        for (int p = 0; p < face_vids[f].size(); ++p)
          face.vertex_ids.push_back(polyh_cell->vertex_ids[face_vids[f][p]]);

        polyh_cell->faces.push_back(face);
      }
      break;
    }
    // We dont need cell vertex id remapping but we need to impose our own
    // face orientation
    case CellType::TETRAHEDRON:
    {
      std::vector<std::vector<uint64_t>> face_vids = {{0, 2, 1}, {0, 1, 3}, {0, 3, 2}, {3, 1, 2}};
      for (int f = 0; f < 4; ++f)
      {
        LightWeightFace face;

        face.vertex_ids.reserve(3);
        for (int p = 0; p < 3; ++p)
          face.vertex_ids.push_back(polyh_cell->vertex_ids[face_vids[f][p]]);

        polyh_cell->faces.push_back(face);
      }
      break;
    }
    default:
    {
      polyh_cell->faces.reserve(num_cfaces);
      for (int f = 0; f < num_cfaces; ++f)
      {
        LightWeightFace face;
        auto vtk_face = vtk_cell->GetFace(f);
        auto num_face_points = vtk_face->GetNumberOfPoints();

        face.vertex_ids.reserve(num_face_points);
        auto face_point_ids = vtk_face->GetPointIds();
        for (int p = 0; p < num_face_points; ++p)
        {
          uint64_t point_id = face_point_ids->GetId(p);
          face.vertex_ids.push_back(point_id);
        }

        polyh_cell->faces.push_back(face);
      }
      break;
    } // default
  }   // switch on sub_type

  return polyh_cell;
}

UnpartitionedMesh::LightWeightCell*
UnpartitionedMesh::CreateCellFromVTKPolygon(vtkCell* vtk_cell)
{
  const std::string fname = "UnpartitionedMesh::CreateCellFromVTKPolygon";

  CellType sub_type;
  switch (vtk_cell->GetCellType())
  {
    case VTK_POLYGON:
      sub_type = CellType::POLYGON;
      break;
    case VTK_QUAD:
    case VTK_PIXEL:
      sub_type = CellType::QUADRILATERAL;
      break;
    case VTK_TRIANGLE:
      sub_type = CellType::TRIANGLE;
      break;
    default:
      throw std::logic_error(fname + ": Unsupported 2D cell type encountered.");
  }

  auto poly_cell = new LightWeightCell(CellType::POLYGON, sub_type);

  auto num_cpoints = vtk_cell->GetNumberOfPoints();
  auto num_cfaces = num_cpoints;

  poly_cell->vertex_ids.reserve(num_cpoints);
  auto point_ids = vtk_cell->GetPointIds();
  for (int p = 0; p < num_cpoints; ++p)
  {
    uint64_t point_id = point_ids->GetId(p);
    poly_cell->vertex_ids.push_back(point_id);
  } // for p

  poly_cell->faces.reserve(num_cfaces);
  for (int f = 0; f < num_cfaces; ++f)
  {
    LightWeightFace face;

    auto v0_id = poly_cell->vertex_ids[f];
    auto v1_id = (f < (num_cfaces - 1)) ? poly_cell->vertex_ids[f + 1] : poly_cell->vertex_ids[0];

    face.vertex_ids.reserve(2);
    face.vertex_ids.push_back(v0_id);
    face.vertex_ids.push_back(v1_id);

    poly_cell->faces.push_back(face);
  }

  return poly_cell;
}

UnpartitionedMesh::LightWeightCell*
UnpartitionedMesh::CreateCellFromVTKLine(vtkCell* vtk_cell)
{
  const std::string fname = "UnpartitionedMesh::CreateCellFromVTKPolygon";

  CellType sub_type;
  switch (vtk_cell->GetCellType())
  {
    case VTK_LINE:
      sub_type = CellType::SLAB;
      break;
    default:
      throw std::logic_error(fname + ": Unsupported 1D cell type encountered.");
  }

  auto slab_cell = new LightWeightCell(CellType::SLAB, sub_type);

  auto vtk_line = vtkLine::SafeDownCast(vtk_cell);
  auto num_cpoints = vtk_line->GetNumberOfPoints();
  auto num_cfaces = num_cpoints;

  slab_cell->vertex_ids.reserve(num_cpoints);
  auto point_ids = vtk_line->GetPointIds();
  for (int p = 0; p < num_cpoints; ++p)
  {
    uint64_t point_id = point_ids->GetId(p);
    slab_cell->vertex_ids.push_back(point_id);
  } // for p

  slab_cell->faces.reserve(num_cfaces);
  for (int f = 0; f < num_cfaces; ++f)
  {
    LightWeightFace face;

    auto v_id = slab_cell->vertex_ids[f];

    face.vertex_ids.reserve(1);
    face.vertex_ids.push_back(v_id);

    slab_cell->faces.push_back(face);
  }

  return slab_cell;
}

UnpartitionedMesh::LightWeightCell*
UnpartitionedMesh::CreateCellFromVTKVertex(vtkCell* vtk_cell)
{
  auto point_cell = new LightWeightCell(CellType::GHOST, CellType::POINT);

  auto vtk_vertex = vtkVertex::SafeDownCast(vtk_cell);
  auto num_cpoints = vtk_vertex->GetNumberOfPoints();

  point_cell->vertex_ids.reserve(num_cpoints);
  auto point_ids = vtk_vertex->GetPointIds();
  for (int p = 0; p < num_cpoints; ++p)
  {
    uint64_t point_id = point_ids->GetId(p);
    point_cell->vertex_ids.push_back(point_id);
  } // for p

  return point_cell;
}

void
UnpartitionedMesh::CopyUGridCellsAndPoints(vtkUnstructuredGrid& ugrid,
                                           const double scale,
                                           int dimension_to_copy)
{
  const std::string fname = "UnpartitionedMesh::CopyUGridCellsAndPoints";
  typedef Vector3 Vec3;
  typedef Vec3* Vec3Ptr;
  typedef LightWeightCell* CellPtr;

  const vtkIdType total_cell_count = ugrid.GetNumberOfCells();
  const vtkIdType total_point_count = ugrid.GetNumberOfPoints();

  bool has_cell_gids = ugrid.GetCellData()->GetGlobalIds();
  bool has_pnts_gids = ugrid.GetPointData()->GetGlobalIds();
  bool has_global_ids = has_cell_gids and has_pnts_gids;

  const auto& block_id_array_name = mesh_options_.material_id_fieldname;

  if (not ugrid.GetCellData()->GetArray(block_id_array_name.c_str()))
    throw std::logic_error(fname + ": grid has no \"" + block_id_array_name + "\" array.");

  auto block_id_array =
    vtkIntArray::SafeDownCast(ugrid.GetCellData()->GetArray(block_id_array_name.c_str()));

  ChiLogicalErrorIf(not block_id_array, "Failed to cast BlockID array to vtkInt");

  if (has_global_ids)
  {
    std::vector<CellPtr> cells(total_cell_count, nullptr);
    std::vector<Vec3Ptr> vertices(total_point_count, nullptr);

    auto cell_gids_ptr = ugrid.GetCellData()->GetGlobalIds();
    auto pnts_gids_ptr = ugrid.GetPointData()->GetGlobalIds();

    auto cell_gids = vtkIdTypeArray::SafeDownCast(cell_gids_ptr);
    auto pnts_gids = vtkIdTypeArray::SafeDownCast(pnts_gids_ptr);

    //=========================================== Determine id offset
    // We do this because some mesh formats (like ExodusII)
    // are indexed with a 1 base instead of 0
    int cid_offset = 0, pid_offset = 0;
    {
      vtkIdType min_cid = total_point_count; // Minimum cell-id
      vtkIdType min_pid = total_point_count; // Minimum point-id

      for (vtkIdType c = 0; c < total_cell_count; ++c)
        min_cid = std::min(min_cid, cell_gids->GetValue(c));

      for (vtkIdType p = 0; p < total_point_count; ++p)
        min_pid = std::min(min_pid, pnts_gids->GetValue(p));

      cid_offset -= static_cast<int>(min_cid);
      pid_offset -= static_cast<int>(min_pid);
    } // build offset

    //=========================================== Build node map
    std::vector<vtkIdType> node_map(total_point_count, 0);
    for (vtkIdType p = 0; p < total_point_count; ++p)
      node_map[p] = pnts_gids->GetValue(p) + pid_offset;

    //=========================================== Load cells
    for (vtkIdType c = 0; c < total_cell_count; ++c)
    {
      auto vtk_cell = ugrid.GetCell(static_cast<vtkIdType>(c));
      auto vtk_celldim = vtk_cell->GetCellDimension();
      const vtkIdType cell_gid = cell_gids->GetValue(c) + cid_offset;

      if (vtk_celldim != dimension_to_copy) continue;

      CellPtr raw_cell;
      if (vtk_celldim == 3) raw_cell = CreateCellFromVTKPolyhedron(vtk_cell);
      else if (vtk_celldim == 2)
        raw_cell = CreateCellFromVTKPolygon(vtk_cell);
      else if (vtk_celldim == 1)
        raw_cell = CreateCellFromVTKLine(vtk_cell);
      else if (vtk_celldim == 0)
        raw_cell = CreateCellFromVTKVertex(vtk_cell);
      else
        throw std::logic_error(fname + ": Unsupported cell dimension ." +
                               std::to_string(vtk_celldim));

      // Map the cell vertex-ids
      for (uint64_t& vid : raw_cell->vertex_ids)
        vid = node_map[vid];

      // Map face vertex-ids
      for (auto& face : raw_cell->faces)
        for (uint64_t& vid : face.vertex_ids)
          vid = node_map[vid];

      raw_cell->material_id = block_id_array->GetValue(c);

      cells[cell_gid] = raw_cell;
    } // for cell c

    //=========================================== Load points
    for (vtkIdType p = 0; p < total_point_count; ++p)
    {
      auto point = ugrid.GetPoint(static_cast<vtkIdType>(p));
      const vtkIdType point_gid = pnts_gids->GetValue(p) + pid_offset;

      auto vertex = new Vec3(point[0], point[1], point[2]);

      *vertex *= scale;

      vertices.at(point_gid) = vertex;
    } // for point p

    //=========================================== Check all cells assigned
    for (vtkIdType c = 0; c < total_cell_count; ++c)
      if (cells[c] == nullptr) throw std::logic_error(fname + ": Cell pointer not assigned ");

    //=========================================== Check all points assigned
    for (vtkIdType p = 0; p < total_point_count; ++p)
      if (vertices[p] == nullptr) throw std::logic_error(fname + ": Vertex pointer not assigned");

    raw_cells_ = cells;
    vertices_.reserve(total_point_count);
    for (auto& vertex_ptr : vertices)
      vertices_.push_back(*vertex_ptr);
  } // If global-ids available
  else
  {
    //======================================== Push cells
    for (vtkIdType c = 0; c < total_cell_count; ++c)
    {
      auto vtk_cell = ugrid.GetCell(static_cast<vtkIdType>(c));
      auto vtk_celldim = vtk_cell->GetCellDimension();

      if (vtk_celldim != dimension_to_copy) continue;

      if (vtk_celldim == 3) raw_cells_.push_back(CreateCellFromVTKPolyhedron(vtk_cell));
      else if (vtk_celldim == 2)
        raw_cells_.push_back(CreateCellFromVTKPolygon(vtk_cell));
      else if (vtk_celldim == 1)
        raw_cells_.push_back(CreateCellFromVTKLine(vtk_cell));
      else if (vtk_celldim == 0)
        raw_cells_.push_back(CreateCellFromVTKVertex(vtk_cell));
      else
        throw std::logic_error(fname + ": Unsupported cell dimension.");

      raw_cells_.back()->material_id = block_id_array->GetValue(c);
    } // for c

    //======================================== Push points
    for (size_t p = 0; p < total_point_count; ++p)
    {
      auto point = ugrid.GetPoint(static_cast<vtkIdType>(p));

      Vec3 vertex(point[0], point[1], point[2]);

      vertex *= scale;

      vertices_.emplace_back(point[0], point[1], point[2]);
    }
  } // if no global-ids

  //================================================== Determine bound box
  for (size_t p = 0; p < total_point_count; ++p)
  {
    const auto& vertex = vertices_[p];
    if (not bound_box_)
      bound_box_ = std::shared_ptr<BoundBox>(
        new BoundBox{vertex.x, vertex.x, vertex.y, vertex.y, vertex.z, vertex.z});

    bound_box_->xmin = std::min(bound_box_->xmin, vertex.x);
    bound_box_->xmax = std::max(bound_box_->xmax, vertex.x);
    bound_box_->ymin = std::min(bound_box_->ymin, vertex.y);
    bound_box_->ymax = std::max(bound_box_->ymax, vertex.y);
    bound_box_->zmin = std::min(bound_box_->zmin, vertex.z);
    bound_box_->zmax = std::max(bound_box_->zmax, vertex.z);
  }

  Chi::log.Log() << fname + ": Done";
}

void
UnpartitionedMesh::SetMaterialIDsFromList(const std::vector<int>& material_ids)
{
  const size_t total_cell_count = raw_cells_.size();

  for (size_t c = 0; c < total_cell_count; ++c)
    raw_cells_[c]->material_id = material_ids[c];
}

void
UnpartitionedMesh::SetBoundaryIDsFromBlocks(std::vector<vtkUGridPtrAndName>& bndry_grid_blocks)
{
  const double EPSILON = 1.0e-12;
  //======================================== Build boundary faces
  std::vector<LightWeightFace*> bndry_faces;
  for (auto& cell_ptr : raw_cells_)
    for (auto& face : cell_ptr->faces)
      if (not face.has_neighbor) bndry_faces.push_back(&face);

  Chi::log.Log() << "Number of boundary faces: " << bndry_faces.size();

  //======================================== Build boundary vertex ids
  std::set<uint64_t> bndry_vids_set;
  for (const auto& face_ptr : bndry_faces)
    for (const auto vid : face_ptr->vertex_ids)
      bndry_vids_set.insert(vid);

  //======================================== Process each boundary block
  uint64_t bndry_id = 0;
  for (const auto& ugrid_name : bndry_grid_blocks)
  {
    auto ugrid = ugrid_name.first;

    mesh_options_.boundary_id_map[bndry_id] = ugrid_name.second;

    //================================= Build vertex map
    bool mapping_failed = false;
    std::vector<size_t> vertex_map(ugrid->GetNumberOfPoints(), 0);
    for (size_t p = 0; p < ugrid->GetNumberOfPoints(); ++p)
    {
      Vector3 point;
      ugrid->GetPoint(static_cast<vtkIdType>(p), &point.x);

      bool map_found = false;
      for (const auto vid : bndry_vids_set)
        if ((point - vertices_[vid]).NormSquare() < EPSILON)
        {
          vertex_map[p] = vid;
          map_found = true;
          break;
        }

      if (not map_found)
      {
        Chi::log.Log0Warning() << "UnpartitionedMesh::"
                                  "SetBoundaryIDsFromBlocks: Failed to map a vertex. " +
                                    point.PrintStr() + " for boundary " + ugrid_name.second +
                                    " therefore the boundary assignment was skipped.";
        mapping_failed = true;
        break;
      }
    } // for point in boundary block

    if (mapping_failed) continue;

    //================================= Build vertex subscriptions
    std::map<uint64_t, std::set<size_t>> vertex_face_subs;
    for (size_t f = 0; f < bndry_faces.size(); ++f)
      for (const auto vid : bndry_faces[f]->vertex_ids)
        vertex_face_subs[vid].insert(f);

    //================================= Process each cell in bndry block
    size_t num_faces_boundarified = 0;
    auto& bndry_block = ugrid_name.first;
    size_t num_bndry_block_cells = bndry_block->GetNumberOfCells();
    for (size_t bc = 0; bc < num_bndry_block_cells; ++bc)
    {
      auto bndry_cell = bndry_block->GetCell(static_cast<vtkIdType>(bc));

      //========================== Build list of face candidates
      //                           and vertex set
      std::set<size_t> face_ids_short_list;
      std::set<uint64_t> bndry_cell_id_set;
      size_t num_points = bndry_cell->GetNumberOfPoints();
      for (size_t p = 0; p < num_points; ++p)
      {
        auto point_id = bndry_cell->GetPointId(static_cast<int>(p));
        bndry_cell_id_set.insert(vertex_map[point_id]);
        for (size_t face_id : vertex_face_subs[vertex_map[point_id]])
          face_ids_short_list.insert(face_id);
      } // for point p

      for (size_t face_id : face_ids_short_list)
      {
        auto& face = bndry_faces[face_id];
        const auto& face_vids = face->vertex_ids;
        std::set<uint64_t> face_id_set(face_vids.begin(), face_vids.end());

        if (face_id_set == bndry_cell_id_set)
        {
          face->neighbor = bndry_id;
          ++num_faces_boundarified;
        }
      } // for face_id
    }   // for boundary cell bc

    Chi::log.Log() << "UnpartitionedMesh: assigned " << num_faces_boundarified << " to boundary id "
                   << bndry_id << " with name " << ugrid_name.second;

    ++bndry_id;
  } // for boundary_block
}

void
UnpartitionedMesh::ReadFromVTU(const UnpartitionedMesh::Options& options)
{
  Chi::log.Log() << "Reading VTU file: " << options.file_name << ".";

  //======================================== Attempt to open file
  std::ifstream file;
  file.open(options.file_name);
  if (!file.is_open()) throw ErrorReadingFile(ReadFromVTU);
  file.close();

  //======================================== Read the file
  mesh_options_ = options;
  auto reader = vtkSmartPointer<vtkXMLUnstructuredGridReader>::New();
  reader->SetFileName(options.file_name.c_str());

  if (not reader->CanReadFile(options.file_name.c_str()))
    throw std::logic_error("Unable to read file-type with this routine");
  reader->UpdateInformation();
  reader->Update();

  //======================================== Get all the grid blocks
  // For vtu files this is very simple. The
  // output of the reader is an UnstructuredGrid.
  auto ugrid_main = vtkUGridPtr(reader->GetOutput());
  std::vector<vtkUGridPtrAndName> grid_blocks = {{ugrid_main, ""}};

  //======================================== Get the main + bndry blocks
  const int max_dimension = FindHighestDimension(grid_blocks);
  Chi::log.Log0Verbose1() << "Maximum dimension : " << max_dimension << "\n";
  std::vector<vtkUGridPtrAndName> domain_grid_blocks =
    GetBlocksOfDesiredDimension(grid_blocks, max_dimension);
  std::vector<vtkUGridPtrAndName> bndry_grid_blocks =
    GetBlocksOfDesiredDimension(grid_blocks, max_dimension - 1);

  //======================================== Process blocks
  auto ugrid = ConsolidateGridBlocks(domain_grid_blocks, mesh_options_.material_id_fieldname);

  //======================================== Copy Data
  CopyUGridCellsAndPoints(*ugrid, options.scale, max_dimension);

  //======================================== Set material ids
  const auto material_ids =
    BuildCellMaterialIDsFromField(ugrid, options.material_id_fieldname, options.file_name);
  SetMaterialIDsFromList(material_ids);

  //======================================== Always do this
  MeshAttributes dimension = NONE;
  switch (max_dimension)
  {
    case 1:
      dimension = DIMENSION_1;
      break;
    case 2:
      dimension = DIMENSION_2;
      break;
    case 3:
      dimension = DIMENSION_3;
      break;
    default:
      break;
  }

  attributes_ = dimension | UNSTRUCTURED;

  ComputeCentroidsAndCheckQuality();
  BuildMeshConnectivity();

  Chi::log.Log() << "Done reading VTU file: " << options.file_name << ".";
}

void
UnpartitionedMesh::ReadFromPVTU(const UnpartitionedMesh::Options& options)
{
  Chi::log.Log() << "Reading PVTU file: " << options.file_name << ".";

  //======================================== Attempt to open file
  std::ifstream file;
  file.open(options.file_name);
  if (!file.is_open()) throw ErrorReadingFile(ReadFromVTU);
  file.close();

  //======================================== Read the file
  mesh_options_ = options;
  auto reader = vtkSmartPointer<vtkXMLPUnstructuredGridReader>::New();
  reader->SetFileName(options.file_name.c_str());

  if (not reader->CanReadFile(options.file_name.c_str()))
    throw std::logic_error("Unable to read file-type with this routine");
  reader->UpdateInformation();
  reader->Update();

  //======================================== Get all the grid blocks
  // For vtu files this is very simple. The
  // output of the reader is an UnstructuredGrid.
  auto ugrid_main = vtkUGridPtr(reader->GetOutput());
  std::vector<vtkUGridPtrAndName> grid_blocks = {{ugrid_main, ""}};

  //======================================== Get the main + bndry blocks
  const int max_dimension = FindHighestDimension(grid_blocks);
  Chi::log.Log0Verbose1() << "Maximum dimension : " << max_dimension << "\n";
  std::vector<vtkUGridPtrAndName> domain_grid_blocks =
    GetBlocksOfDesiredDimension(grid_blocks, max_dimension);
  std::vector<vtkUGridPtrAndName> bndry_grid_blocks =
    GetBlocksOfDesiredDimension(grid_blocks, max_dimension - 1);

  //======================================== Process blocks
  auto ugrid = ConsolidateGridBlocks(domain_grid_blocks);

  //======================================== Copy Data
  CopyUGridCellsAndPoints(*ugrid, options.scale, max_dimension);

  //======================================== Set material ids
  const auto material_ids =
    BuildCellMaterialIDsFromField(ugrid, options.material_id_fieldname, options.file_name);
  SetMaterialIDsFromList(material_ids);

  //======================================== Always do this
  MeshAttributes dimension = NONE;
  switch (max_dimension)
  {
    case 1:
      dimension = DIMENSION_1;
      break;
    case 2:
      dimension = DIMENSION_2;
      break;
    case 3:
      dimension = DIMENSION_3;
      break;
    default:
      break;
  }

  attributes_ = dimension | UNSTRUCTURED;

  ComputeCentroidsAndCheckQuality();
  BuildMeshConnectivity();

  Chi::log.Log() << "Done reading PVTU file: " << options.file_name << ".";
}

void
UnpartitionedMesh::ReadFromEnsightGold(const UnpartitionedMesh::Options& options)
{
  Chi::log.Log() << "Reading Ensight-Gold file: " << options.file_name << ".";

  //======================================== Attempt to open file
  std::ifstream file;
  file.open(options.file_name);
  if (!file.is_open()) throw ErrorReadingFile(ReadFromEnsightGold);
  file.close();

  //======================================== Read the file
  mesh_options_ = options;
  auto reader = vtkSmartPointer<vtkEnSightGoldBinaryReader>::New();
  reader->SetCaseFileName(options.file_name.c_str());

  if (not reader->CanReadFile(options.file_name.c_str()))
    throw std::logic_error("Unable to read file-type with this routine");
  reader->UpdateInformation();
  reader->Update();

  //======================================== Get all the grid blocks
  auto multiblock = reader->GetOutput();

  std::vector<vtkUGridPtrAndName> grid_blocks;
  auto iter_a = multiblock->NewIterator();
  iter_a->GoToFirstItem();
  while (not iter_a->IsDoneWithTraversal())
  {
    auto block_a = iter_a->GetCurrentDataObject();

    const std::string block_name =
      chi::StringTrim(iter_a->GetCurrentMetaData()->Get(vtkCompositeDataSet::NAME()));

    if (block_a->GetDataObjectType() == VTK_UNSTRUCTURED_GRID)
    {
      grid_blocks.emplace_back(
        vtkUnstructuredGrid::SafeDownCast(block_a),
        chi::StringTrim(iter_a->GetCurrentMetaData()->Get(vtkCompositeDataSet::NAME())));

      Chi::log.Log() << "Reading block " << block_name
                     << " Number of cells: " << grid_blocks.back().first->GetNumberOfCells()
                     << " Number of points: " << grid_blocks.back().first->GetNumberOfPoints();
    }

    iter_a->GoToNextItem();
  }

  //======================================== Get the main + bndry blocks
  const int max_dimension = FindHighestDimension(grid_blocks);
  Chi::log.Log0Verbose1() << "Maximum dimension : " << max_dimension << "\n";
  std::vector<vtkUGridPtrAndName> domain_grid_blocks =
    GetBlocksOfDesiredDimension(grid_blocks, max_dimension);
  std::vector<vtkUGridPtrAndName> bndry_grid_blocks =
    GetBlocksOfDesiredDimension(grid_blocks, max_dimension - 1);

  //======================================== Process blocks
  SetBlockIDArrays(domain_grid_blocks);
  auto ugrid = ConsolidateGridBlocks(domain_grid_blocks);

  //======================================== Copy Data
  // Material-IDs will get set form block-id arrays
  CopyUGridCellsAndPoints(*ugrid, options.scale, max_dimension);

  //======================================== Always do this
  MeshAttributes dimension = NONE;
  switch (max_dimension)
  {
    case 1:
      dimension = DIMENSION_1;
      break;
    case 2:
      dimension = DIMENSION_2;
      break;
    case 3:
      dimension = DIMENSION_3;
      break;
    default:
      break;
  }

  attributes_ = dimension | UNSTRUCTURED;

  ComputeCentroidsAndCheckQuality();
  BuildMeshConnectivity();

  //======================================== Set boundary ids
  SetBoundaryIDsFromBlocks(bndry_grid_blocks);

  Chi::log.Log() << "Done reading Ensight-Gold file: " << options.file_name << ".";
}

void
UnpartitionedMesh::ReadFromWavefrontOBJ(const Options& options)
{
  const std::string fname = "UnpartitionedMesh::ReadFromWavefrontOBJ";

  //======================================================= Opening the file
  std::ifstream file;
  file.open(options.file_name);
  if (!file.is_open())
  {
    Chi::log.LogAllError() << "Failed to open file: " << options.file_name << " in call "
                           << "to ImportFromOBJFile \n";
    Chi::Exit(EXIT_FAILURE);
  }

  Chi::mpi.Barrier();
  Chi::log.Log() << "Making Unpartitioned mesh from wavefront file " << options.file_name;

  typedef std::pair<uint64_t, uint64_t> Edge;
  struct BlockData
  {
    std::string name;
    std::vector<LightWeightCell*> cells;
    std::vector<Edge> edges;
  };

  std::vector<BlockData> block_data;
  std::vector<Vertex> file_vertices;

  //======================================================= Reading every line
  std::string file_line;
  std::string delimiter = " ";
  int material_id = -1;
  while (std::getline(file, file_line))
  {
    //================================================ Get the first word
    size_t beg_of_word = file_line.find_first_not_of(delimiter);
    size_t end_of_word = file_line.find(delimiter, beg_of_word - beg_of_word);
    std::string first_word = file_line.substr(beg_of_word, end_of_word);
    std::string sub_word;

    if (first_word == "o")
    {
      beg_of_word = file_line.find_first_not_of(delimiter, end_of_word);
      end_of_word = file_line.find(delimiter, beg_of_word);
      sub_word = file_line.substr(beg_of_word, end_of_word - beg_of_word);

      std::string block_name = sub_word;
      block_data.push_back({block_name, {}});
    }

    if (first_word == "usemtl")
    {
      Chi::log.Log0Verbose1() << "New material at cell count: " << block_data.back().cells.size();
      ++material_id;
    }

    //================================================ Keyword "v" for Vertex
    if (first_word == "v")
    {
      Vertex newVertex;
      for (int k = 1; k <= 3; k++)
      {
        //================================== Extract sub word
        beg_of_word = file_line.find_first_not_of(delimiter, end_of_word);
        end_of_word = file_line.find(delimiter, beg_of_word);
        sub_word = file_line.substr(beg_of_word, end_of_word - beg_of_word);

        //================================== Convert word to number
        try
        {
          double numValue = std::stod(sub_word);

          if (k == 1) newVertex.x = numValue;
          else if (k == 2)
            newVertex.y = numValue;
          else if (k == 3)
            newVertex.z = numValue;
        }

        //================================== Catch conversion error
        catch (const std::invalid_argument& ia)
        {
          Chi::log.Log0Warning() << "Failed to convert vertex in line " << file_line << std::endl;
        }

        //================================== Stop word extraction on line end
        if (end_of_word == std::string::npos) { break; }
      }
      file_vertices.push_back(newVertex);
    } // if (first_word == "v")

    //===================================================== Keyword "f" for face
    if (first_word == "f")
    {
      size_t number_of_verts = std::count(file_line.begin(), file_line.end(), '/') / 2;

      CellType sub_type = CellType::POLYGON;
      if (number_of_verts == 3) sub_type = CellType::TRIANGLE;
      else if (number_of_verts == 4)
        sub_type = CellType::QUADRILATERAL;

      auto cell = new LightWeightCell(CellType::POLYGON, sub_type);
      cell->material_id = material_id;

      // Populate vertex-ids
      for (size_t k = 1; k <= number_of_verts; k++)
      {
        //================================== Extract sub word
        beg_of_word = file_line.find_first_not_of(delimiter, end_of_word);
        end_of_word = file_line.find(delimiter, beg_of_word);
        sub_word = file_line.substr(beg_of_word, end_of_word - beg_of_word);

        //============================= Extract locations of hiphens
        size_t first_dash = sub_word.find('/');
        size_t last_dash = sub_word.find_last_of('/');

        //============================= Extract the words ass. vertex and normal
        std::string vert_word = sub_word.substr(0, first_dash - 0);
        std::string norm_word = sub_word.substr(last_dash + 1, sub_word.length() - last_dash - 1);

        //============================= Convert word to number (Vertex)
        try
        {
          int numValue = std::stoi(vert_word);
          cell->vertex_ids.push_back(numValue - 1);
        }
        catch (const std::invalid_argument& ia)
        {
          Chi::log.Log0Warning() << "Failed converting work to number in line " << file_line
                                 << std::endl;
        }

        //============================= Stop word extraction on line end
        if (end_of_word == std::string::npos) { break; }
      }

      // Build faces
      const size_t num_verts = cell->vertex_ids.size();
      for (uint64_t v = 0; v < num_verts; ++v)
      {
        LightWeightFace face;

        face.vertex_ids.resize(2);
        face.vertex_ids[0] = cell->vertex_ids[v];
        face.vertex_ids[1] = (v < (num_verts - 1)) ? cell->vertex_ids[v + 1] : cell->vertex_ids[0];

        cell->faces.push_back(std::move(face));
      } // for v

      if (block_data.empty())
        throw std::logic_error(fname + ": Could not add cell to block-data. "
                                       "This normally indicates that the file does not have the "
                                       "\"o Object Name\" entry.");

      block_data.back().cells.push_back(cell);
    } // if (first_word == "f")

    //===================================================== Keyword "l" for edge
    if (first_word == "l")
    {
      Edge edge;
      for (int k = 1; k <= 2; ++k)
      {
        //================================== Extract sub word
        beg_of_word = file_line.find_first_not_of(delimiter, end_of_word);
        end_of_word = file_line.find(delimiter, beg_of_word);
        sub_word = file_line.substr(beg_of_word, end_of_word - beg_of_word);

        //================================== Convert word to number
        try
        {
          int vertex_id = std::stoi(sub_word);
          if (k == 1) edge.first = vertex_id - 1;
          if (k == 2) edge.second = vertex_id - 1;
        }

        //================================== Catch conversion error
        catch (const std::invalid_argument& ia)
        {
          Chi::log.Log0Warning() << "Failed to text to integer in line " << file_line << std::endl;
        }
      } // for k

      if (block_data.empty())
        throw std::logic_error(fname + ": Could not add edge to block-data. "
                                       "This normally indicates that the file does not have the "
                                       "\"o Object Name\" entry.");

      block_data.back().edges.push_back(edge);
    } // if (first_word == "l")
  }
  file.close();
  Chi::log.Log0Verbose0() << "Max material id: " << material_id;

  //======================================================= Error checks
  for (const auto& block : block_data)
    for (const auto& cell : block.cells)
      for (const auto vid : cell->vertex_ids)
      {
        ChiLogicalErrorIf(
          vid >= file_vertices.size(),
          "Cell vertex id " + std::to_string(vid) +
            " not in list of vertices read (size=" + std::to_string(file_vertices.size()) + ").");
      }

  //======================================================= Filter blocks
  std::vector<size_t> bndry_block_ids;
  size_t num_cell_blocks = 0;
  size_t main_block_id = 0;
  for (size_t block_id = 0; block_id < block_data.size(); ++block_id)
  {
    if (not block_data[block_id].edges.empty()) bndry_block_ids.push_back(block_id);
    if (not block_data[block_id].cells.empty())
    {
      ++num_cell_blocks;
      main_block_id = block_id;
    }
  } // for block_id

  if (num_cell_blocks != 1)
    throw std::logic_error(fname +
                           ": More than one cell-block has been read "
                           "from the file. Only a single face-containing object is supported. "
                           "If you exported this mesh from blender, be sure to export "
                           "\"selection only\"");

  //======================================================= Process blocks
  std::vector<Vertex> cell_vertices;
  {
    // Initial map is straight
    std::vector<size_t> vertex_map;
    vertex_map.reserve(file_vertices.size());
    for (size_t m = 0; m < file_vertices.size(); ++m)
      vertex_map.push_back(m);

    // Build set of cell vertices
    std::set<size_t> cell_vertex_id_set;
    for (const auto& cell_ptr : block_data.at(main_block_id).cells)
      for (size_t vid : cell_ptr->vertex_ids)
        cell_vertex_id_set.insert(vid);

    // Make cell_vertices and edit map
    {
      size_t new_id = 0;
      for (size_t vid : cell_vertex_id_set)
      {
        cell_vertices.push_back(file_vertices[vid]);
        vertex_map[vid] = new_id;
        ++new_id;
      }
    }

    // Build set of bndry vertices
    std::set<size_t> bndry_vertex_id_set;
    for (size_t block_id : bndry_block_ids)
      for (const auto& edge : block_data[block_id].edges)
      {
        bndry_vertex_id_set.insert(edge.first);
        bndry_vertex_id_set.insert(edge.second);
      }

    // Find a match for each boundary vertex and
    // place it in the map
    for (size_t bvid : bndry_vertex_id_set)
    {
      const auto& bndry_vertex = file_vertices[bvid];

      bool match_found = false;
      for (size_t cvid = 0; cvid < cell_vertices.size(); ++cvid)
      {
        const auto& cell_vertex = cell_vertices[cvid];

        if ((bndry_vertex - cell_vertex).NormSquare() < 1.0e-12)
        {
          vertex_map[bvid] = cvid;
          match_found = true;
          break;
        }
      } // for cvid
      if (not match_found)
        throw std::logic_error(fname +
                               ": Failed to map a boundary vertex to"
                               "any cell vertex. Check that the edges are conformal with the "
                               "object containing the faces.");
    } // for bvid

    // Change cell and face vertex ids to cell vertex ids
    // using vertex map
    for (auto& cell_ptr : block_data.at(main_block_id).cells)
    {
      for (uint64_t& vid : cell_ptr->vertex_ids)
        vid = vertex_map[vid];

      for (auto& face : cell_ptr->faces)
        for (uint64_t& vid : face.vertex_ids)
          vid = vertex_map[vid];
    }

    // Change edge vertex ids to cell vertex ids using
    // the vertex map
    for (size_t block_id : bndry_block_ids)
      for (auto& edge : block_data[block_id].edges)
      {
        edge.first = static_cast<int>(vertex_map[edge.first]);
        edge.second = static_cast<int>(vertex_map[edge.second]);
      }
  }
  this->vertices_ = cell_vertices;
  this->raw_cells_ = block_data[main_block_id].cells;

  //======================================================= Always do this
  attributes_ = DIMENSION_2 | UNSTRUCTURED;

  ComputeCentroidsAndCheckQuality();
  BuildMeshConnectivity();

  //======================================================= Set boundary ids
  if (bndry_block_ids.empty()) { mesh_options_.boundary_id_map[0] = "Default Boundary"; }
  else
  {
    std::vector<LightWeightFace*> bndry_faces;
    for (auto& cell_ptr : raw_cells_)
      for (auto& face : cell_ptr->faces)
        if (not face.has_neighbor) bndry_faces.push_back(&face);

    size_t bndry_id = 0;
    for (size_t bid : bndry_block_ids)
    {
      const auto& bndry_edges = block_data[bid].edges;

      size_t num_faces_boundarified = 0;
      for (const auto& edge : bndry_edges)
      {
        std::set<size_t> edge_vert_id_set({edge.first, edge.second});

        for (auto& face_ptr : bndry_faces)
        {
          const auto& vert_ids = face_ptr->vertex_ids;
          std::set<size_t> face_vert_id_set(vert_ids.begin(), vert_ids.end());

          if (face_vert_id_set == edge_vert_id_set)
          {
            face_ptr->neighbor = bndry_id;
            ++num_faces_boundarified;
            break;
          }
        } // for face
      }   // for edge

      Chi::log.Log() << "UnpartitionedMesh: assigned " << num_faces_boundarified
                     << " faces to boundary id " << bndry_id << " with name "
                     << block_data[bid].name;

      mesh_options_.boundary_id_map[bndry_id] = block_data[bid].name;

      ++bndry_id;
    } // for boundary block
  }
}

void
UnpartitionedMesh::ReadFromMsh(const Options& options)
{
  const std::string fname = "UnpartitionedMesh::ReadFromMsh";

  //===================================================== Opening the file
  std::ifstream file;
  file.open(options.file_name);
  if (!file.is_open())
  {
    Chi::log.LogAllError() << "Failed to open file: " << options.file_name << " in call "
                           << "to ReadFromMsh \n";
    Chi::Exit(EXIT_FAILURE);
  }

  Chi::log.Log() << "Making Unpartitioned mesh from msh format file " << options.file_name;
  Chi::mpi.Barrier();

  //===================================================== Declarations
  std::string file_line;
  std::istringstream iss;
  const std::string node_section_name = "$Nodes";
  const std::string elements_section_name = "$Elements";
  const std::string format_section_name = "$MeshFormat";

  //=================================================== Check the format of this
  // input
  // Role file forward until "$MeshFormat" line is encountered.
  while (std::getline(file, file_line))
    if (format_section_name == file_line) break;

  std::getline(file, file_line);
  iss = std::istringstream(file_line);
  double format;
  if (!(iss >> format)) throw std::logic_error(fname + ": Failed to read the file format.");
  else if (format != 2.2)
    throw std::logic_error(fname + ": Currently, only msh format 2.2 is supported.");

  //=================================================== Find section with node
  //                                                    information and then
  //                                                    read the nodes
  file.seekg(0);
  while (std::getline(file, file_line))
    if (node_section_name == file_line) break;

  std::getline(file, file_line);
  iss = std::istringstream(file_line);
  int num_nodes;
  if (!(iss >> num_nodes))
    throw std::logic_error(fname + ": Failed while trying to read "
                                   "the number of nodes.");

  vertices_.clear();
  vertices_.resize(num_nodes);

  for (int n = 0; n < num_nodes; n++)
  {
    std::getline(file, file_line);
    iss = std::istringstream(file_line);

    int vert_index;
    if (!(iss >> vert_index)) throw std::logic_error(fname + ": Failed to read vertex index.");

    if (!(iss >> vertices_[vert_index - 1].x >> vertices_[vert_index - 1].y >>
          vertices_[vert_index - 1].z))
      throw std::logic_error(fname + ": Failed while reading the vertex "
                                     "coordinates.");
  }

  //================================================== Define utility lambdas
  /**Lambda for reading nodes.*/
  auto ReadNodes = [&iss, &fname](int num_nodes)
  {
    std::vector<int> raw_nodes(num_nodes, 0);
    for (int i = 0; i < num_nodes; ++i)
      if (!(iss >> raw_nodes[i]))
        throw std::logic_error(fname + ": Failed when reading element "
                                       "node index.");

    std::vector<uint64_t> nodes(num_nodes, 0);
    for (int i = 0; i < num_nodes; ++i)
      if ((raw_nodes[i] - 1) >= 0) nodes[i] = raw_nodes[i] - 1;
    return nodes;
  };

  /**Lamda for checking if an element is 1D.*/
  auto IsElementType1D = [](int element_type)
  {
    if (element_type == 1) return true;

    return false;
  };

  /**Lamda for checking if an element is 2D.*/
  auto IsElementType2D = [](int element_type)
  {
    if (element_type == 2 or element_type == 3) return true;

    return false;
  };

  /**Lamda for checking if an element is 2D.*/
  auto IsElementType3D = [](int element_type)
  {
    if (element_type >= 4 and element_type <= 7) return true;

    return false;
  };

  /**Lambda for checking supported elements.*/
  auto IsElementSupported = [](int element_type)
  {
    if (element_type >= 1 and element_type <= 7) return true;

    return false;
  };

  /**Lambda giving the cell subtype, given the MSH cell type.*/
  auto CellTypeFromMSHTypeID = [](int element_type)
  {
    CellType cell_type = CellType::GHOST;

    if (element_type == 1) cell_type = CellType::SLAB;
    else if (element_type == 2)
      cell_type = CellType::TRIANGLE;
    else if (element_type == 3)
      cell_type = CellType::QUADRILATERAL;
    else if (element_type == 4)
      cell_type = CellType::TETRAHEDRON;
    else if (element_type == 5)
      cell_type = CellType::HEXAHEDRON;
    else if (element_type == 6 or // Prism
             element_type == 7)
      cell_type = CellType::POLYHEDRON; // Pyramid

    return cell_type;
  };

  //================================================== Determine mesh type 2D/3D
  // Only 2D and 3D meshes are supported. If the mesh
  // is 1D then no elements will be read but the state
  // would still be safe.
  // This section will run through all the elements
  // looking for a 3D element. It will not process
  // any elements.
  bool mesh_is_2D_assumption = true;
  file.seekg(0);
  while (std::getline(file, file_line))
    if (elements_section_name == file_line) break;

  std::getline(file, file_line);
  iss = std::istringstream(file_line);
  int num_elems;
  if (!(iss >> num_elems)) throw std::logic_error(fname + ": Failed to read number of elements.");

  for (int n = 0; n < num_elems; n++)
  {
    int elem_type, num_tags, physical_reg, tag, element_index;

    std::getline(file, file_line);
    iss = std::istringstream(file_line);

    if (!(iss >> element_index >> elem_type >> num_tags))
      throw std::logic_error(fname + ": Failed while reading element index, "
                                     "element type, and number of tags.");

    if (!(iss >> physical_reg))
      throw std::logic_error(fname + ": Failed while reading physical region.");

    for (int i = 1; i < num_tags; i++)
      if (!(iss >> tag)) throw std::logic_error(fname + ": Failed when reading tags.");

    if (IsElementType3D(elem_type))
    {
      mesh_is_2D_assumption = false;
      Chi::log.Log() << "Mesh identified as 3D.";
      break; // have the answer now leave loop
    }

    if (elem_type == 15) // skip point type element
      continue;

    if (not IsElementSupported(elem_type))
      throw std::logic_error(fname + ": Unsupported element encountered.");
  } // for n

  //================================================== Return to the element
  //                                                   listing section
  // Now we will actually read the elements.
  file.seekg(0);
  while (std::getline(file, file_line))
    if (elements_section_name == file_line) break;

  std::getline(file, file_line);
  iss = std::istringstream(file_line);
  if (!(iss >> num_elems)) throw std::logic_error(fname + ": Failed to read number of elements.");

  for (int n = 0; n < num_elems; n++)
  {
    int elem_type, num_tags, physical_reg, tag, element_index;

    std::getline(file, file_line);
    iss = std::istringstream(file_line);

    if (!(iss >> element_index >> elem_type >> num_tags))
      throw std::logic_error(fname + ": Failed while reading element index, "
                                     "element type, and number of tags.");

    if (!(iss >> physical_reg))
      throw std::logic_error(fname + ": Failed while reading physical region.");

    for (int i = 1; i < num_tags; i++)
      if (!(iss >> tag)) throw std::logic_error(fname + ": Failed when reading tags.");

    if (elem_type == 15) // skip point type elements
      continue;

    if (not IsElementSupported(elem_type))
      throw std::logic_error(fname + ": Unsupported element encountered.");

    Chi::log.Log0Verbose2() << "Reading element: " << file_line << " type: " << elem_type;

    int num_cell_nodes;
    if (elem_type == 1) num_cell_nodes = 2;
    else if (elem_type == 2) // 3-node triangle
      num_cell_nodes = 3;
    else if (elem_type == 3 or elem_type == 4) // 4-node quadrangle or tet
      num_cell_nodes = 4;
    else if (elem_type == 5) // 8-node hexahedron
      num_cell_nodes = 8;
    else
      continue;

    //====================================== Make the cell on either the volume
    //                                       or the boundary
    LightWeightCell* raw_cell = nullptr;
    if (mesh_is_2D_assumption)
    {
      if (IsElementType1D(elem_type))
      {
        raw_cell = new LightWeightCell(CellType::SLAB, CellType::SLAB);
        raw_boundary_cells_.push_back(raw_cell);
        Chi::log.Log0Verbose2() << "Added to raw_boundary_cells.";
      }
      else if (IsElementType2D(elem_type))
      {
        raw_cell = new LightWeightCell(CellType::POLYGON, CellTypeFromMSHTypeID(elem_type));
        raw_cells_.push_back(raw_cell);
        Chi::log.Log0Verbose2() << "Added to raw_cells.";
      }
    }
    else
    {
      if (IsElementType2D(elem_type))
      {
        raw_cell = new LightWeightCell(CellType::POLYGON, CellTypeFromMSHTypeID(elem_type));
        raw_boundary_cells_.push_back(raw_cell);
        Chi::log.Log0Verbose2() << "Added to raw_boundary_cells.";
      }
      else if (IsElementType3D(elem_type))
      {
        raw_cell = new LightWeightCell(CellType::POLYHEDRON, CellTypeFromMSHTypeID(elem_type));
        raw_cells_.push_back(raw_cell);
        Chi::log.Log0Verbose2() << "Added to raw_cells.";
      }
    }

    if (raw_cell == nullptr) continue;

    auto& cell = *raw_cell;
    cell.material_id = physical_reg;
    cell.vertex_ids = ReadNodes(num_cell_nodes);

    //====================================== Populate faces
    if (elem_type == 1) // 2-node edge
    {
      LightWeightFace face0;
      LightWeightFace face1;

      face0.vertex_ids = {cell.vertex_ids.at(0)};
      face1.vertex_ids = {cell.vertex_ids.at(1)};

      cell.faces.push_back(face0);
      cell.faces.push_back(face1);
    }
    else if (elem_type == 2 or elem_type == 3) // 3-node triangle or 4-node quadrangle
    {
      size_t num_verts = cell.vertex_ids.size();
      for (size_t e = 0; e < num_verts; e++)
      {
        size_t ep1 = (e < (num_verts - 1)) ? e + 1 : 0;
        LightWeightFace face;

        face.vertex_ids = {cell.vertex_ids[e], cell.vertex_ids[ep1]};

        cell.faces.push_back(std::move(face));
      }                      // for e
    }                        // if 2D elements
    else if (elem_type == 4) // 4-node tetrahedron
    {
      auto& v = cell.vertex_ids;
      std::vector<LightWeightFace> lw_faces(4);
      lw_faces[0].vertex_ids = {v[0], v[2], v[1]}; // base-face
      lw_faces[1].vertex_ids = {v[0], v[3], v[2]};
      lw_faces[2].vertex_ids = {v[3], v[1], v[2]};
      lw_faces[3].vertex_ids = {v[3], v[0], v[1]};

      for (auto& lw_face : lw_faces)
        cell.faces.push_back(lw_face);
    }
    else if (elem_type == 5) // 8-node hexahedron
    {
      auto& v = cell.vertex_ids;
      std::vector<LightWeightFace> lw_faces(6);
      lw_faces[0].vertex_ids = {v[5], v[1], v[2], v[6]}; // East face
      lw_faces[1].vertex_ids = {v[0], v[4], v[7], v[3]}; // West face
      lw_faces[2].vertex_ids = {v[0], v[3], v[2], v[1]}; // North face
      lw_faces[3].vertex_ids = {v[4], v[5], v[6], v[7]}; // South face
      lw_faces[4].vertex_ids = {v[2], v[3], v[7], v[6]}; // Top face
      lw_faces[5].vertex_ids = {v[0], v[1], v[5], v[4]}; // Bottom face

      for (auto& lw_face : lw_faces)
        cell.faces.push_back(lw_face);
    }
    else
      throw std::runtime_error(fname + ": Unsupported cell type");

  } // for elements

  file.close();

  //======================================== Remap material-ids
  std::set<int> material_ids_set_as_read;
  std::map<int, int> material_mapping;

  for (auto& cell : raw_cells_)
    material_ids_set_as_read.insert(cell->material_id);

  std::set<int> boundary_ids_set_as_read;
  std::map<int, int> boundary_mapping;

  for (auto& cell : raw_boundary_cells_)
    boundary_ids_set_as_read.insert(cell->material_id);

  {
    int m = 0;
    for (const auto& mat_id : material_ids_set_as_read)
      material_mapping.insert(std::make_pair(mat_id, m++));

    int b = 0;
    for (const auto& bndry_id : boundary_ids_set_as_read)
      boundary_mapping.insert(std::make_pair(bndry_id, b++));
  }

  for (auto& cell : raw_cells_)
    cell->material_id = material_mapping[cell->material_id];

  for (auto& cell : raw_boundary_cells_)
    cell->material_id = boundary_mapping[cell->material_id];

  //======================================== Always do this
  MeshAttributes dimension = DIMENSION_2;
  if (not mesh_is_2D_assumption) dimension = DIMENSION_3;

  attributes_ = dimension | UNSTRUCTURED;

  ComputeCentroidsAndCheckQuality();
  BuildMeshConnectivity();

  Chi::log.Log() << "Done processing " << options.file_name << ".\n"
                 << "Number of nodes read: " << vertices_.size() << "\n"
                 << "Number of cells read: " << raw_cells_.size();
}

/**Reads an Exodus unstructured mesh.*/
void
UnpartitionedMesh::ReadFromExodus(const UnpartitionedMesh::Options& options)
{
  Chi::log.Log() << "Reading Exodus file: " << options.file_name << ".";

  //======================================== Attempt to open file
  std::ifstream file;
  file.open(options.file_name);
  if (!file.is_open()) throw ErrorReadingFile(ReadFromExodus);
  file.close();

  //======================================== Read the file
  mesh_options_ = options;
  auto reader = vtkSmartPointer<vtkExodusIIReader>::New();
  reader->SetFileName(options.file_name.c_str());

  if (not reader->CanReadFile(options.file_name.c_str()))
    throw std::logic_error("Unable to read file-type with this routine");

  reader->UpdateInformation();
  // Exodus ships boundary-ids via SideSets and NodeSets. This allows
  // it to be read from the file. Here we have to enable the reader
  // to process this because it does not do it by default.
  reader->SetAllArrayStatus(reader->NODE_SET, 1);
  reader->SetAllArrayStatus(reader->NODE_SET_CONN, 1);
  reader->SetAllArrayStatus(reader->SIDE_SET, 1);
  reader->SetAllArrayStatus(reader->SIDE_SET_CONN, 1);

  // The exodusII file format ships blocks of elements
  // together with their points/vertices as self-contained (localized)
  // unstructured meshes. To relate these localized blocks to the original
  // mesh, where are the blocks formed a whole, we need to know the mapping
  // from block-local ids to the original ids. This information can
  // be derived from the GlobalNodeID arrays loaded onto point-data and
  // cell-data. Again, this information is not read by default so we have to
  // turn this on.
  reader->SetGenerateGlobalNodeIdArray(true);
  reader->SetGenerateGlobalElementIdArray(true);
  reader->Update();

  //======================================== Get all the grid blocks
  // This part was quite difficult. I eventually found how to do
  // this from this post:
  // https://public.kitware.com/pipermail/vtkusers/2010-December/064921.html
  // It indicated that all exodus formats are read with a single
  // top level vtkMultiBlockDataSet (Level 1). Each of the blocks in
  // level 2 is also of type vtkMultiBlockDataSet. The level 2 block that has
  // the actual elements is also split into blocks but these, level 3,
  // blocks each contain a structure castable to vtkUnstructuredGrid.
  auto multiblock = reader->GetOutput();

  std::vector<vtkUGridPtrAndName> grid_blocks;
  auto iter_a = multiblock->NewIterator();
  iter_a->GoToFirstItem();
  while (not iter_a->IsDoneWithTraversal())
  {
    auto block_a = iter_a->GetCurrentDataObject();

    const std::string block_name =
      chi::StringTrim(iter_a->GetCurrentMetaData()->Get(vtkCompositeDataSet::NAME()));

    if (block_a->GetDataObjectType() == VTK_UNSTRUCTURED_GRID)
    {
      grid_blocks.emplace_back(vtkUnstructuredGrid::SafeDownCast(block_a), block_name);

      Chi::log.Log() << "Reading block " << block_name
                     << " Number of cells: " << grid_blocks.back().first->GetNumberOfCells()
                     << " Number of points: " << grid_blocks.back().first->GetNumberOfPoints();
    }

    iter_a->GoToNextItem();
  }

  //======================================== Get the main + bndry blocks
  const int max_dimension = FindHighestDimension(grid_blocks);
  Chi::log.Log0Verbose1() << "Maximum dimension : " << max_dimension << "\n";
  std::vector<vtkUGridPtrAndName> domain_grid_blocks =
    GetBlocksOfDesiredDimension(grid_blocks, max_dimension);
  std::vector<vtkUGridPtrAndName> bndry_grid_blocks =
    GetBlocksOfDesiredDimension(grid_blocks, max_dimension - 1);

  //======================================== Process blocks
  SetBlockIDArrays(domain_grid_blocks);
  auto ugrid = ConsolidateGridBlocks(domain_grid_blocks);

  //======================================== Copy Data
  // Material-IDs will get set form block-id arrays
  CopyUGridCellsAndPoints(*ugrid, options.scale, max_dimension);

  //======================================== Always do this
  MeshAttributes dimension = NONE;
  switch (max_dimension)
  {
    case 1:
      dimension = DIMENSION_1;
      break;
    case 2:
      dimension = DIMENSION_2;
      break;
    case 3:
      dimension = DIMENSION_3;
      break;
    default:
      break;
  }

  attributes_ = dimension | UNSTRUCTURED;

  ComputeCentroidsAndCheckQuality();
  BuildMeshConnectivity();

  //======================================== Set boundary ids
  SetBoundaryIDsFromBlocks(bndry_grid_blocks);

  Chi::log.Log() << "Done reading Exodus file: " << options.file_name << ".";
}

/**Makes a cell from proxy information and pushes the cell to the mesh.*/
void
UnpartitionedMesh::PushProxyCell(const std::string& type_str,
                                 const std::string& sub_type_str,
                                 int cell_num_faces,
                                 int cell_material_id,
                                 const std::vector<std::vector<uint64_t>>& proxy_faces)
{
  const std::string fname = __FUNCTION__;
  CellType type;
  if (type_str == "SLAB") type = CellType::SLAB;
  else if (type_str == "POLYGON")
    type = CellType::POLYGON;
  else if (type_str == "POLYHEDRON")
    type = CellType::POLYHEDRON;
  else
    throw std::logic_error(fname + ": Unsupported cell primary type.");

  CellType sub_type;
  if (sub_type_str == "SLAB") sub_type = CellType::SLAB;
  else if (sub_type_str == "POLYGON")
    sub_type = CellType::POLYGON;
  else if (sub_type_str == "TRIANGLE")
    sub_type = CellType::TRIANGLE;
  else if (sub_type_str == "QUADRILATERAL")
    sub_type = CellType::QUADRILATERAL;
  else if (sub_type_str == "POLYHEDRON")
    sub_type = CellType::POLYHEDRON;
  else if (sub_type_str == "TETRAHEDRON")
    sub_type = CellType::TETRAHEDRON;
  else if (sub_type_str == "HEXAHEDRON")
    sub_type = CellType::HEXAHEDRON;
  else
    throw std::logic_error(fname + ": Unsupported cell secondary type.");

  auto cell = new LightWeightCell(type, sub_type);

  cell->material_id = cell_material_id;

  // Filter cell-vertex-ids from faces
  std::set<uint64_t> cell_vertex_id_set;
  for (auto& proxy_face : proxy_faces)
    for (uint64_t fvid : proxy_face)
      cell_vertex_id_set.insert(fvid);

  // Assign cell-vertex-ids
  cell->vertex_ids.reserve(cell_vertex_id_set.size());
  for (uint64_t cvid : cell_vertex_id_set)
    cell->vertex_ids.push_back(cvid);

  // Assign faces from proxy faces
  cell->faces.reserve(cell_num_faces);
  for (auto& proxy_face : proxy_faces)
    cell->faces.emplace_back(proxy_face);

  raw_cells_.push_back(cell);

  MeshAttributes dimension;
  if (type == CellType::SLAB) dimension = DIMENSION_1;
  if (type == CellType::POLYGON) dimension = DIMENSION_2;
  if (type == CellType::POLYHEDRON) dimension = DIMENSION_3;

  attributes_ = dimension | UNSTRUCTURED;
}

} // namespace chi_mesh
