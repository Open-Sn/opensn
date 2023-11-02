#include "chi_unpartitioned_mesh.h"

#include <vtkUnstructuredGrid.h>
#include <vtkAppendFilter.h>

#include <vtkCellData.h>
#include <vtkPointData.h>

#include "chi_runtime.h"
#include "chi_log.h"

// ###################################################################
/**Copies the vtk data structures to the current object's internal
 * data.*/
void
chi_mesh::UnpartitionedMesh::CopyUGridCellsAndPoints(vtkUnstructuredGrid& ugrid,
                                                     const double scale,
                                                     int dimension_to_copy)
{
  const std::string fname = "chi_mesh::UnpartitionedMesh::CopyUGridCellsAndPoints";
  typedef chi_mesh::Vector3 Vec3;
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

// ###################################################################
/**Set material-ids from list.*/
void
chi_mesh::UnpartitionedMesh::SetMaterialIDsFromList(const std::vector<int>& material_ids)
{
  const size_t total_cell_count = raw_cells_.size();

  for (size_t c = 0; c < total_cell_count; ++c)
    raw_cells_[c]->material_id = material_ids[c];
}

// ###################################################################
/**Set boundary-ids from boundary grid_blocks.*/
void
chi_mesh::UnpartitionedMesh::SetBoundaryIDsFromBlocks(
  std::vector<vtkUGridPtrAndName>& bndry_grid_blocks)
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
      chi_mesh::Vector3 point;
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
        Chi::log.Log0Warning() << "chi_mesh::UnpartitionedMesh::"
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
