// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "framework/io/vtk_io.h"
#include "framework/runtime.h"
#include "framework/logging/log.h"
#include "framework/utils/utils.h"
#include <vtkCell.h>
#include <vtkPolygon.h>
#include <vtkLine.h>
#include <vtkVertex.h>
#include <vtkUnstructuredGrid.h>
#include <vtkCellData.h>
#include <vtkPointData.h>
#include <vtkSmartPointer.h>
#include <vtkXMLUnstructuredGridReader.h>
#include <vtkExodusIIReader.h>
#include <vtkInformation.h>
#include <vtkXMLPUnstructuredGridReader.h>
#include <vtkEnSightGoldBinaryReader.h>
#include <vtkMultiBlockDataSet.h>
#include <fstream>

namespace opensn
{

namespace
{

UnpartitionedMesh::LightWeightCell*
CreateCellFromVTKPolyhedron(vtkCell* vtk_cell)
{
  const std::string fname = "CreateCellFromVTKPolyhedron";

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
  auto polyh_cell = new UnpartitionedMesh::LightWeightCell(CellType::POLYHEDRON, sub_type);

  auto num_cpoints = vtk_cell->GetNumberOfPoints();
  auto num_cfaces = vtk_cell->GetNumberOfFaces();

  polyh_cell->vertex_ids.reserve(num_cpoints);
  auto point_ids = vtk_cell->GetPointIds();
  for (int p = 0; p < num_cpoints; ++p)
  {
    uint64_t point_id = point_ids->GetId(p);
    polyh_cell->vertex_ids.push_back(point_id);
  }

  switch (sub_type)
  {
    // The cell vertex ids in VTK is the same as in OpenSn, so we don't
    // need to remap the vertices. We do however need to remap the faces.
    case CellType::HEXAHEDRON:
    {
      std::vector<std::vector<uint64_t>> face_vids = {
        {1, 2, 6, 5}, {3, 0, 4, 7}, {2, 3, 7, 6}, {0, 1, 5, 4}, {4, 5, 6, 7}, {3, 2, 1, 0}};
      for (int f = 0; f < 6; ++f)
      {
        UnpartitionedMesh::LightWeightFace face;

        face.vertex_ids.reserve(4);
        for (int p = 0; p < 4; ++p)
          face.vertex_ids.push_back(polyh_cell->vertex_ids[face_vids[f][p]]);

        polyh_cell->faces.push_back(face);
      }
      break;
    }
    // For wedges, we need to remap cell vertices and faces.
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
        UnpartitionedMesh::LightWeightFace face;

        face.vertex_ids.reserve(4);
        for (int p = 0; p < face_vids[f].size(); ++p)
          face.vertex_ids.push_back(polyh_cell->vertex_ids[face_vids[f][p]]);

        polyh_cell->faces.push_back(face);
      }
      break;
    }
    // We don't need cell vertex id remapping, but we need to impose our own face orientation
    case CellType::TETRAHEDRON:
    {
      std::vector<std::vector<uint64_t>> face_vids = {{0, 2, 1}, {0, 1, 3}, {0, 3, 2}, {3, 1, 2}};
      for (int f = 0; f < 4; ++f)
      {
        UnpartitionedMesh::LightWeightFace face;

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
        UnpartitionedMesh::LightWeightFace face;
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
    }
  }

  return polyh_cell;
}

UnpartitionedMesh::LightWeightCell*
CreateCellFromVTKPolygon(vtkCell* vtk_cell)
{
  const std::string fname = "CreateCellFromVTKPolygon";

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

  auto poly_cell = new UnpartitionedMesh::LightWeightCell(CellType::POLYGON, sub_type);

  auto num_cpoints = vtk_cell->GetNumberOfPoints();
  auto num_cfaces = num_cpoints;

  poly_cell->vertex_ids.reserve(num_cpoints);
  auto point_ids = vtk_cell->GetPointIds();
  for (int p = 0; p < num_cpoints; ++p)
  {
    uint64_t point_id = point_ids->GetId(p);
    poly_cell->vertex_ids.push_back(point_id);
  }

  poly_cell->faces.reserve(num_cfaces);
  for (int f = 0; f < num_cfaces; ++f)
  {
    UnpartitionedMesh::LightWeightFace face;

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
CreateCellFromVTKLine(vtkCell* vtk_cell)
{
  const std::string fname = "CreateCellFromVTKPolygon";

  CellType sub_type;
  switch (vtk_cell->GetCellType())
  {
    case VTK_LINE:
      sub_type = CellType::SLAB;
      break;
    default:
      throw std::logic_error(fname + ": Unsupported 1D cell type encountered.");
  }

  auto slab_cell = new UnpartitionedMesh::LightWeightCell(CellType::SLAB, sub_type);

  auto vtk_line = vtkLine::SafeDownCast(vtk_cell);
  auto num_cpoints = vtk_line->GetNumberOfPoints();
  auto num_cfaces = num_cpoints;

  slab_cell->vertex_ids.reserve(num_cpoints);
  auto point_ids = vtk_line->GetPointIds();
  for (int p = 0; p < num_cpoints; ++p)
  {
    uint64_t point_id = point_ids->GetId(p);
    slab_cell->vertex_ids.push_back(point_id);
  }

  slab_cell->faces.reserve(num_cfaces);
  for (int f = 0; f < num_cfaces; ++f)
  {
    UnpartitionedMesh::LightWeightFace face;

    auto v_id = slab_cell->vertex_ids[f];

    face.vertex_ids.reserve(1);
    face.vertex_ids.push_back(v_id);

    slab_cell->faces.push_back(face);
  }

  return slab_cell;
}

UnpartitionedMesh::LightWeightCell*
CreateCellFromVTKVertex(vtkCell* vtk_cell)
{
  auto point_cell = new UnpartitionedMesh::LightWeightCell(CellType::GHOST, CellType::POINT);

  auto vtk_vertex = vtkVertex::SafeDownCast(vtk_cell);
  auto num_cpoints = vtk_vertex->GetNumberOfPoints();

  point_cell->vertex_ids.reserve(num_cpoints);
  auto point_ids = vtk_vertex->GetPointIds();
  for (int p = 0; p < num_cpoints; ++p)
  {
    uint64_t point_id = point_ids->GetId(p);
    point_cell->vertex_ids.push_back(point_id);
  }

  return point_cell;
}

void
CopyUGridCellsAndPoints(std::shared_ptr<UnpartitionedMesh> mesh,
                        vtkUnstructuredGrid& ugrid,
                        const double scale,
                        int dimension_to_copy,
                        const std::string& block_id_array_name)
{
  const std::string fname = "CopyUGridCellsAndPoints";

  const vtkIdType total_cell_count = ugrid.GetNumberOfCells();
  const vtkIdType total_point_count = ugrid.GetNumberOfPoints();

  bool has_cell_gids = ugrid.GetCellData()->GetGlobalIds();
  bool has_pnts_gids = ugrid.GetPointData()->GetGlobalIds();
  bool has_global_ids = has_cell_gids and has_pnts_gids;

  if (not ugrid.GetCellData()->GetArray(block_id_array_name.c_str()))
    throw std::logic_error(fname + ": grid has no \"" + block_id_array_name + "\" array.");

  auto block_id_array =
    vtkIntArray::SafeDownCast(ugrid.GetCellData()->GetArray(block_id_array_name.c_str()));

  OpenSnLogicalErrorIf(not block_id_array, "Failed to cast BlockID array to vtkInt");

  if (has_global_ids)
  {
    std::vector<UnpartitionedMesh::LightWeightCell*> cells(total_cell_count, nullptr);
    std::vector<Vector3*> vertices(total_point_count, nullptr);

    auto cell_gids_ptr = ugrid.GetCellData()->GetGlobalIds();
    auto pnts_gids_ptr = ugrid.GetPointData()->GetGlobalIds();

    auto cell_gids = vtkIdTypeArray::SafeDownCast(cell_gids_ptr);
    auto pnts_gids = vtkIdTypeArray::SafeDownCast(pnts_gids_ptr);

    // Determine id offset
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

    // Build node map
    std::vector<vtkIdType> node_map(total_point_count, 0);
    for (vtkIdType p = 0; p < total_point_count; ++p)
      node_map[p] = pnts_gids->GetValue(p) + pid_offset;

    // Load cells
    for (vtkIdType c = 0; c < total_cell_count; ++c)
    {
      auto vtk_cell = ugrid.GetCell(static_cast<vtkIdType>(c));
      auto vtk_celldim = vtk_cell->GetCellDimension();
      const vtkIdType cell_gid = cell_gids->GetValue(c) + cid_offset;

      if (vtk_celldim != dimension_to_copy)
        continue;

      UnpartitionedMesh::LightWeightCell* raw_cell;
      if (vtk_celldim == 3)
        raw_cell = CreateCellFromVTKPolyhedron(vtk_cell);
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

    // Load points
    for (vtkIdType p = 0; p < total_point_count; ++p)
    {
      auto point = ugrid.GetPoint(static_cast<vtkIdType>(p));
      const vtkIdType point_gid = pnts_gids->GetValue(p) + pid_offset;

      auto vertex = new Vector3(point[0], point[1], point[2]);

      *vertex *= scale;

      vertices.at(point_gid) = vertex;
    } // for point p

    // Check all cells assigned
    for (vtkIdType c = 0; c < total_cell_count; ++c)
      if (cells[c] == nullptr)
        throw std::logic_error(fname + ": Cell pointer not assigned ");

    // Check all points assigned
    for (vtkIdType p = 0; p < total_point_count; ++p)
      if (vertices[p] == nullptr)
        throw std::logic_error(fname + ": Vertex pointer not assigned");

    mesh->RawCells() = cells;
    mesh->Vertices().reserve(total_point_count);
    for (auto& vertex_ptr : vertices)
      mesh->Vertices().push_back(*vertex_ptr);
  } // If global-ids available
  else
  {
    auto& raw_cells = mesh->RawCells();
    // Push cells
    for (vtkIdType c = 0; c < total_cell_count; ++c)
    {
      auto vtk_cell = ugrid.GetCell(static_cast<vtkIdType>(c));
      auto vtk_celldim = vtk_cell->GetCellDimension();

      if (vtk_celldim != dimension_to_copy)
        continue;

      if (vtk_celldim == 3)
        raw_cells.push_back(CreateCellFromVTKPolyhedron(vtk_cell));
      else if (vtk_celldim == 2)
        raw_cells.push_back(CreateCellFromVTKPolygon(vtk_cell));
      else if (vtk_celldim == 1)
        raw_cells.push_back(CreateCellFromVTKLine(vtk_cell));
      else if (vtk_celldim == 0)
        raw_cells.push_back(CreateCellFromVTKVertex(vtk_cell));
      else
        throw std::logic_error(fname + ": Unsupported cell dimension.");

      raw_cells.back()->material_id = block_id_array->GetValue(c);
    }

    // Push points
    for (size_t p = 0; p < total_point_count; ++p)
    {
      auto point = ugrid.GetPoint(static_cast<vtkIdType>(p));

      Vector3 vertex(point[0], point[1], point[2]);

      vertex *= scale;

      mesh->Vertices().emplace_back(point[0], point[1], point[2]);
    }
  } // if no global-ids

  mesh->ComputeBoundingBox();

  log.Log() << fname + ": Done";
}

void
SetMaterialIDsFromList(std::shared_ptr<UnpartitionedMesh> mesh,
                       const std::vector<int>& material_ids)
{
  auto& raw_cells = mesh->RawCells();
  const size_t total_cell_count = raw_cells.size();
  for (size_t c = 0; c < total_cell_count; ++c)
    raw_cells[c]->material_id = material_ids[c];
}

void
SetBoundaryIDsFromBlocks(std::shared_ptr<UnpartitionedMesh> mesh,
                         std::vector<vtkUGridPtrAndName>& bndry_grid_blocks)
{
  const double EPSILON = 1.0e-12;
  auto& raw_cells = mesh->RawCells();
  // Build boundary faces
  std::vector<UnpartitionedMesh::LightWeightFace*> bndry_faces;
  for (auto& cell_ptr : raw_cells)
    for (auto& face : cell_ptr->faces)
      if (not face.has_neighbor)
        bndry_faces.push_back(&face);

  log.Log() << "Number of boundary faces: " << bndry_faces.size();

  // Build boundary vertex ids
  std::set<uint64_t> bndry_vids_set;
  for (const auto& face_ptr : bndry_faces)
    for (const auto vid : face_ptr->vertex_ids)
      bndry_vids_set.insert(vid);

  // Process each boundary block
  uint64_t bndry_id = 0;
  for (const auto& ugrid_name : bndry_grid_blocks)
  {
    auto ugrid = ugrid_name.first;

    mesh->AddBoundary(bndry_id, ugrid_name.second);

    // Build vertex map
    bool mapping_failed = false;
    std::vector<size_t> vertex_map(ugrid->GetNumberOfPoints(), 0);
    for (size_t p = 0; p < ugrid->GetNumberOfPoints(); ++p)
    {
      Vector3 point;
      ugrid->GetPoint(static_cast<vtkIdType>(p), &point.x);

      bool map_found = false;
      for (const auto vid : bndry_vids_set)
        if ((point - mesh->Vertices()[vid]).NormSquare() < EPSILON)
        {
          vertex_map[p] = vid;
          map_found = true;
          break;
        }

      if (not map_found)
      {
        log.Log0Warning() << "UnpartitionedMesh::"
                             "SetBoundaryIDsFromBlocks: Failed to map a vertex. " +
                               point.PrintStr() + " for boundary " + ugrid_name.second +
                               " therefore the boundary assignment was skipped.";
        mapping_failed = true;
        break;
      }
    } // for point in boundary block

    if (mapping_failed)
      continue;

    // Build vertex subscriptions
    std::map<uint64_t, std::set<size_t>> vertex_face_subs;
    for (size_t f = 0; f < bndry_faces.size(); ++f)
      for (const auto vid : bndry_faces[f]->vertex_ids)
        vertex_face_subs[vid].insert(f);

    // Process each cell in bndry block
    size_t num_faces_boundarified = 0;
    auto& bndry_block = ugrid_name.first;
    size_t num_bndry_block_cells = bndry_block->GetNumberOfCells();
    for (size_t bc = 0; bc < num_bndry_block_cells; ++bc)
    {
      auto bndry_cell = bndry_block->GetCell(static_cast<vtkIdType>(bc));

      // Build list of face candidates and vertex set
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

    log.Log() << "UnpartitionedMesh: assigned " << num_faces_boundarified << " to boundary id "
              << bndry_id << " with name " << ugrid_name.second;

    ++bndry_id;
  } // for boundary_block
}

} // namespace

std::shared_ptr<UnpartitionedMesh>
VtkIO::FromExodus(const UnpartitionedMesh::Options& options)
{
  log.Log() << "Reading Exodus file: " << options.file_name << ".";

  std::shared_ptr<UnpartitionedMesh> mesh = std::make_shared<UnpartitionedMesh>();

  // Read the file
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

  // Get all the grid blocks
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
      StringTrim(iter_a->GetCurrentMetaData()->Get(vtkCompositeDataSet::NAME()));

    if (block_a->GetDataObjectType() == VTK_UNSTRUCTURED_GRID)
    {
      grid_blocks.emplace_back(vtkUnstructuredGrid::SafeDownCast(block_a), block_name);

      log.Log() << "Reading block " << block_name
                << " Number of cells: " << grid_blocks.back().first->GetNumberOfCells()
                << " Number of points: " << grid_blocks.back().first->GetNumberOfPoints();
    }

    iter_a->GoToNextItem();
  }

  // Get the main + bndry blocks
  const int max_dimension = FindHighestDimension(grid_blocks);
  log.Log0Verbose1() << "Maximum dimension : " << max_dimension << "\n";
  std::vector<vtkUGridPtrAndName> domain_grid_blocks =
    GetBlocksOfDesiredDimension(grid_blocks, max_dimension);
  std::vector<vtkUGridPtrAndName> bndry_grid_blocks =
    GetBlocksOfDesiredDimension(grid_blocks, max_dimension - 1);

  // Process blocks
  SetBlockIDArrays(domain_grid_blocks);
  auto ugrid = ConsolidateGridBlocks(domain_grid_blocks);

  // Copy Data
  // Material-IDs will get set form block-id arrays
  CopyUGridCellsAndPoints(
    mesh, *ugrid, options.scale, max_dimension, options.material_id_fieldname);

  // Always do this
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

  mesh->Attributes() = dimension | UNSTRUCTURED;

  mesh->ComputeCentroids();
  mesh->CheckQuality();
  mesh->BuildMeshConnectivity();

  // Set boundary ids
  SetBoundaryIDsFromBlocks(mesh, bndry_grid_blocks);

  log.Log() << "Done reading Exodus file: " << options.file_name << ".";

  return mesh;
}

std::shared_ptr<UnpartitionedMesh>
VtkIO::FromVTU(const UnpartitionedMesh::Options& options)
{
  log.Log() << "Reading VTU file: " << options.file_name << ".";

  std::shared_ptr<UnpartitionedMesh> mesh = std::make_shared<UnpartitionedMesh>();

  // Read the file
  auto reader = vtkSmartPointer<vtkXMLUnstructuredGridReader>::New();
  reader->SetFileName(options.file_name.c_str());

  if (not reader->CanReadFile(options.file_name.c_str()))
    throw std::logic_error("Unable to read file-type with this routine");
  reader->UpdateInformation();
  reader->Update();

  // Get all the grid blocks
  // For vtu files this is very simple. The
  // output of the reader is an UnstructuredGrid.
  auto ugrid_main = vtkUGridPtr(reader->GetOutput());
  std::vector<vtkUGridPtrAndName> grid_blocks = {{ugrid_main, ""}};

  // Get the main + bndry blocks
  const int max_dimension = FindHighestDimension(grid_blocks);
  log.Log0Verbose1() << "Maximum dimension : " << max_dimension << "\n";
  std::vector<vtkUGridPtrAndName> domain_grid_blocks =
    GetBlocksOfDesiredDimension(grid_blocks, max_dimension);
  std::vector<vtkUGridPtrAndName> bndry_grid_blocks =
    GetBlocksOfDesiredDimension(grid_blocks, max_dimension - 1);

  // Process blocks
  auto ugrid = ConsolidateGridBlocks(domain_grid_blocks, options.material_id_fieldname);

  // Copy Data
  CopyUGridCellsAndPoints(
    mesh, *ugrid, options.scale, max_dimension, options.material_id_fieldname);

  // Set material ids
  const auto material_ids =
    BuildCellMaterialIDsFromField(ugrid, options.material_id_fieldname, options.file_name);
  SetMaterialIDsFromList(mesh, material_ids);

  // Always do this
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

  mesh->Attributes() = dimension | UNSTRUCTURED;

  mesh->ComputeCentroids();
  mesh->CheckQuality();
  mesh->BuildMeshConnectivity();

  log.Log() << "Done reading VTU file: " << options.file_name << ".";

  return mesh;
}

std::shared_ptr<UnpartitionedMesh>
VtkIO::FromPVTU(const UnpartitionedMesh::Options& options)
{
  log.Log() << "Reading PVTU file: " << options.file_name << ".";

  std::shared_ptr<UnpartitionedMesh> mesh = std::make_shared<UnpartitionedMesh>();

  // Read the file
  auto reader = vtkSmartPointer<vtkXMLPUnstructuredGridReader>::New();
  reader->SetFileName(options.file_name.c_str());

  if (not reader->CanReadFile(options.file_name.c_str()))
    throw std::logic_error("Unable to read file-type with this routine");
  reader->UpdateInformation();
  reader->Update();

  // Get all the grid blocks
  // For vtu files this is very simple. The
  // output of the reader is an UnstructuredGrid.
  auto ugrid_main = vtkUGridPtr(reader->GetOutput());
  std::vector<vtkUGridPtrAndName> grid_blocks = {{ugrid_main, ""}};

  // Get the main + bndry blocks
  const int max_dimension = FindHighestDimension(grid_blocks);
  log.Log0Verbose1() << "Maximum dimension : " << max_dimension << "\n";
  std::vector<vtkUGridPtrAndName> domain_grid_blocks =
    GetBlocksOfDesiredDimension(grid_blocks, max_dimension);
  std::vector<vtkUGridPtrAndName> bndry_grid_blocks =
    GetBlocksOfDesiredDimension(grid_blocks, max_dimension - 1);

  // Process blocks
  auto ugrid = ConsolidateGridBlocks(domain_grid_blocks);

  // Copy Data
  CopyUGridCellsAndPoints(
    mesh, *ugrid, options.scale, max_dimension, options.material_id_fieldname);

  // Set material ids
  const auto material_ids =
    BuildCellMaterialIDsFromField(ugrid, options.material_id_fieldname, options.file_name);
  SetMaterialIDsFromList(mesh, material_ids);

  // Always do this
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

  mesh->Attributes() = dimension | UNSTRUCTURED;

  mesh->ComputeCentroids();
  mesh->CheckQuality();
  mesh->BuildMeshConnectivity();

  log.Log() << "Done reading PVTU file: " << options.file_name << ".";

  return mesh;
}

std::shared_ptr<UnpartitionedMesh>
VtkIO::FromEnsightGold(const UnpartitionedMesh::Options& options)
{
  log.Log() << "Reading Ensight-Gold file: " << options.file_name << ".";

  std::shared_ptr<UnpartitionedMesh> mesh = std::make_shared<UnpartitionedMesh>();

  // Read the file
  auto reader = vtkSmartPointer<vtkEnSightGoldBinaryReader>::New();
  reader->SetCaseFileName(options.file_name.c_str());

  if (not reader->CanReadFile(options.file_name.c_str()))
    throw std::logic_error("Unable to read file-type with this routine");
  reader->UpdateInformation();
  reader->Update();

  // Get all the grid blocks
  auto multiblock = reader->GetOutput();

  std::vector<vtkUGridPtrAndName> grid_blocks;
  auto iter_a = multiblock->NewIterator();
  iter_a->GoToFirstItem();
  while (not iter_a->IsDoneWithTraversal())
  {
    auto block_a = iter_a->GetCurrentDataObject();

    const std::string block_name =
      StringTrim(iter_a->GetCurrentMetaData()->Get(vtkCompositeDataSet::NAME()));

    if (block_a->GetDataObjectType() == VTK_UNSTRUCTURED_GRID)
    {
      grid_blocks.emplace_back(
        vtkUnstructuredGrid::SafeDownCast(block_a),
        StringTrim(iter_a->GetCurrentMetaData()->Get(vtkCompositeDataSet::NAME())));

      log.Log() << "Reading block " << block_name
                << " Number of cells: " << grid_blocks.back().first->GetNumberOfCells()
                << " Number of points: " << grid_blocks.back().first->GetNumberOfPoints();
    }

    iter_a->GoToNextItem();
  }

  // Get the main + bndry blocks
  const int max_dimension = FindHighestDimension(grid_blocks);
  log.Log0Verbose1() << "Maximum dimension : " << max_dimension << "\n";
  std::vector<vtkUGridPtrAndName> domain_grid_blocks =
    GetBlocksOfDesiredDimension(grid_blocks, max_dimension);
  std::vector<vtkUGridPtrAndName> bndry_grid_blocks =
    GetBlocksOfDesiredDimension(grid_blocks, max_dimension - 1);

  // Process blocks
  SetBlockIDArrays(domain_grid_blocks);
  auto ugrid = ConsolidateGridBlocks(domain_grid_blocks);

  // Copy Data
  // Material-IDs will get set form block-id arrays
  CopyUGridCellsAndPoints(
    mesh, *ugrid, options.scale, max_dimension, options.material_id_fieldname);

  // Always do this
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

  mesh->Attributes() = dimension | UNSTRUCTURED;

  mesh->ComputeCentroids();
  mesh->CheckQuality();
  mesh->BuildMeshConnectivity();

  // Set boundary ids
  SetBoundaryIDsFromBlocks(mesh, bndry_grid_blocks);

  log.Log() << "Done reading Ensight-Gold file: " << options.file_name << ".";

  return mesh;
}

} // namespace opensn
