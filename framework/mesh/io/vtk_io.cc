// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "framework/mesh/io/mesh_io.h"
#include "framework/mesh/mesh_continuum/mesh_continuum.h"
#include "framework/runtime.h"
#include "framework/logging/log.h"
#include "framework/utils/utils.h"
#include "framework/mesh/mesh_continuum/grid_vtk_utils.h"
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
#include <vtkExodusIIWriter.h>
#include <vtkModelMetadata.h>
#include <fstream>

namespace opensn
{

namespace
{

std::shared_ptr<UnpartitionedMesh::LightWeightCell>
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
  auto polyh_cell =
    std::make_shared<UnpartitionedMesh::LightWeightCell>(CellType::POLYHEDRON, sub_type);

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
    // For wedges, we need to remap cell faces
    case CellType::WEDGE:
    {
      std::vector<std::vector<uint64_t>> face_vids = {
        {0, 1, 4, 3}, {1, 2, 5, 4}, {2, 0, 3, 5}, {3, 4, 5}, {0, 2, 1}};
      for (int f = 0; f < 5; ++f)
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

std::shared_ptr<UnpartitionedMesh::LightWeightCell>
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

  auto poly_cell =
    std::make_shared<UnpartitionedMesh::LightWeightCell>(CellType::POLYGON, sub_type);

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

std::shared_ptr<UnpartitionedMesh::LightWeightCell>
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

  auto slab_cell = std::make_shared<UnpartitionedMesh::LightWeightCell>(CellType::SLAB, sub_type);

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

std::shared_ptr<UnpartitionedMesh::LightWeightCell>
CreateCellFromVTKVertex(vtkCell* vtk_cell)
{
  auto point_cell =
    std::make_shared<UnpartitionedMesh::LightWeightCell>(CellType::GHOST, CellType::POINT);

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
    std::vector<std::shared_ptr<UnpartitionedMesh::LightWeightCell>> cells(total_cell_count);
    std::vector<std::shared_ptr<Vector3>> vertices(total_point_count);

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

      std::shared_ptr<UnpartitionedMesh::LightWeightCell> raw_cell;
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

      auto vertex = std::make_shared<Vector3>(point[0], point[1], point[2]);

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

    mesh->GetRawCells() = cells;
    mesh->GetVertices().reserve(total_point_count);
    for (auto& vertex_ptr : vertices)
      mesh->GetVertices().push_back(*vertex_ptr);
  } // If global-ids available
  else
  {
    auto& raw_cells = mesh->GetRawCells();
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

      mesh->GetVertices().emplace_back(point[0], point[1], point[2]);
    }
  } // if no global-ids

  mesh->ComputeBoundingBox();

  log.Log() << fname + ": Done";
}

void
SetMaterialIDsFromList(std::shared_ptr<UnpartitionedMesh> mesh,
                       const std::vector<int>& material_ids)
{
  auto& raw_cells = mesh->GetRawCells();
  const size_t total_cell_count = raw_cells.size();
  for (size_t c = 0; c < total_cell_count; ++c)
    raw_cells[c]->material_id = material_ids[c];
}

void
SetBoundaryIDsFromBlocks(std::shared_ptr<UnpartitionedMesh> mesh,
                         std::vector<vtkUGridPtrAndName>& bndry_grid_blocks)
{
  const double EPSILON = 1.0e-12;
  auto& raw_cells = mesh->GetRawCells();
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
        if ((point - mesh->GetVertices()[vid]).NormSquare() < EPSILON)
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
MeshIO::FromExodusII(const UnpartitionedMesh::Options& options)
{
  log.Log() << "Reading ExodusII file: " << options.file_name << ".";

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
  mesh->SetDimension(max_dimension);
  mesh->SetType(UNSTRUCTURED);

  mesh->ComputeCentroids();
  mesh->CheckQuality();
  mesh->BuildMeshConnectivity();

  // Set boundary ids
  SetBoundaryIDsFromBlocks(mesh, bndry_grid_blocks);

  log.Log() << "Done reading ExodusII file: " << options.file_name << ".";

  return mesh;
}

std::shared_ptr<UnpartitionedMesh>
MeshIO::FromVTU(const UnpartitionedMesh::Options& options)
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
  mesh->SetDimension(max_dimension);
  mesh->SetType(UNSTRUCTURED);

  mesh->ComputeCentroids();
  mesh->CheckQuality();
  mesh->BuildMeshConnectivity();

  log.Log() << "Done reading VTU file: " << options.file_name << ".";

  return mesh;
}

std::shared_ptr<UnpartitionedMesh>
MeshIO::FromPVTU(const UnpartitionedMesh::Options& options)
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
  mesh->SetDimension(max_dimension);
  mesh->SetType(UNSTRUCTURED);

  mesh->ComputeCentroids();
  mesh->CheckQuality();
  mesh->BuildMeshConnectivity();

  log.Log() << "Done reading PVTU file: " << options.file_name << ".";

  return mesh;
}

std::shared_ptr<UnpartitionedMesh>
MeshIO::FromEnsightGold(const UnpartitionedMesh::Options& options)
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
  mesh->SetDimension(max_dimension);
  mesh->SetType(UNSTRUCTURED);

  mesh->ComputeCentroids();
  mesh->CheckQuality();
  mesh->BuildMeshConnectivity();

  // Set boundary ids
  SetBoundaryIDsFromBlocks(mesh, bndry_grid_blocks);

  log.Log() << "Done reading Ensight-Gold file: " << options.file_name << ".";

  return mesh;
}

void
MeshIO::ToOBJ(const std::shared_ptr<MeshContinuum>& grid, const char* file_name, bool per_material)
{
  if (not per_material)
  {
    FILE* of = fopen(file_name, "w");

    if (of == nullptr)
      throw std::logic_error("Could not open file '" + std::string(file_name) + "' for writing.");

    // Develop list of faces and nodes
    std::set<int> nodes_set;
    std::vector<CellFace> faces_to_export;
    for (auto& cell : grid->local_cells)
    {
      if (cell.GetType() == CellType::POLYHEDRON)
      {
        for (auto& face : cell.faces)
        {
          if (not face.has_neighbor)
          {
            faces_to_export.push_back(face);

            for (int vid : face.vertex_ids)
              nodes_set.insert(vid);
          } // if boundary
        }   // for face
      }     // if polyhedron
    }       // for local cell

    // Write header
    fprintf(of, "# Exported mesh file from Extrusion script\n");
    std::string str_file_name(file_name);
    std::string file_base_name = str_file_name.substr(0, str_file_name.find('.'));
    fprintf(of, "o %s\n", file_base_name.c_str());

    // Develop node mapping and write them
    std::vector<int> node_mapping(grid->GetGlobalVertexCount(), -1);

    int node_counter = 0;
    for (auto node : nodes_set)
    {
      node_counter++;
      int node_g_index = node;
      node_mapping[node_g_index] = node_counter;

      Vector3 cur_v = grid->vertices[node_g_index];

      fprintf(of, "v %9.6f %9.6f %9.6f\n", cur_v.x, cur_v.y, cur_v.z);
    }

    // Write face normals
    for (const auto& face : faces_to_export)
    {
      fprintf(of, "vn %.4f %.4f %.4f\n", face.normal.x, face.normal.y, face.normal.z);
    }

    // Write faces
    int normal_counter = 0;
    for (const auto& face : faces_to_export)
    {
      normal_counter++;
      fprintf(of, "f");

      for (auto v_g_index : face.vertex_ids)
        fprintf(of, " %d//%d", node_mapping[v_g_index], normal_counter);

      fprintf(of, "\n");
    }

    fclose(of);

    log.Log() << "Exported Volume mesh to " << str_file_name;
  } // Whole mesh
  // PER MATERIAL
  else
  {
    // Get base name
    std::string str_file_name(file_name);
    std::string file_base_name = str_file_name.substr(0, str_file_name.find('.'));

    if (material_stack.empty())
    {
      log.Log0Warning()
        << "MeshIO::ToOBJ: No mesh will be exported because there are no physics materials present";
    }

    for (int mat = 0; mat < material_stack.size(); ++mat)
    {
      std::string mat_base_name = file_base_name + std::string("_m") + std::to_string(mat);
      std::string mat_file_name = mat_base_name + std::string(".obj");
      FILE* of = fopen(mat_file_name.c_str(), "w");

      if (of == nullptr)
        throw std::logic_error("Could not open file '" + mat_file_name + "' for writing.");

      // Develop list of faces and nodes
      std::set<int> nodes_set;
      std::vector<CellFace> faces_to_export;
      for (const auto& cell : grid->local_cells)
      {
        if (cell.GetType() == CellType::POLYHEDRON)
        {
          if (cell.material_id != mat)
            continue;

          for (const auto& face : cell.faces)
          {
            int adjcell_glob_index = face.neighbor_id;

            if (adjcell_glob_index < 0)
            {
              faces_to_export.push_back(face);

              for (auto vid : face.vertex_ids)
                nodes_set.insert(vid);
            } // if boundary
            else
            {
              auto& adj_cell = grid->cells[adjcell_glob_index];

              if (adj_cell.material_id != mat)
              {
                faces_to_export.push_back(face);

                for (auto vid : face.vertex_ids)
                  nodes_set.insert(vid);
              } // if material mismatch
            }   // if neighbor cell
          }     // for face
        }       // if polyhedron
      }         // for local cell

      // Write header
      fprintf(of, "# Exported mesh file from Extrusion script\n");
      fprintf(of, "o %s\n", mat_base_name.c_str());

      // Develop node mapping and write them
      std::vector<int> node_mapping(grid->GetGlobalVertexCount(), -1);

      int node_counter = 0;
      for (auto node : nodes_set)
      {
        node_counter++;
        int node_g_index = node;
        node_mapping[node_g_index] = node_counter;

        Vector3 cur_v = grid->vertices[node_g_index];

        fprintf(of, "v %9.6f %9.6f %9.6f\n", cur_v.x, cur_v.y, cur_v.z);
      }

      // Write face normals
      for (const auto& face : faces_to_export)
      {
        fprintf(of, "vn %.4f %.4f %.4f\n", face.normal.x, face.normal.y, face.normal.z);
      }

      // Write faces
      int normal_counter = 0;
      for (const auto& face : faces_to_export)
      {
        normal_counter++;
        fprintf(of, "f");

        for (auto v_g_index : face.vertex_ids)
          fprintf(of, " %d//%d", node_mapping[v_g_index], normal_counter);

        fprintf(of, "\n");
      }

      fclose(of);

      log.Log() << "Exported Material Volume mesh to " << mat_file_name;
    } // for mat
  }   // if per material
  opensn::mpi_comm.barrier();
}

void
MeshIO::ToExodusII(const std::shared_ptr<MeshContinuum>& grid,
                   const std::string& file_name,
                   bool write_node_sets,
                   bool write_side_sets)
{
  const std::string fname = "MeshIO::ToExodusII";
  log.Log() << "Exporting mesh to ExodusII file with base " << file_name;

  if (opensn::mpi_comm.size() != 1)
    throw std::logic_error(fname + ": Currently this routine is only allowed in serial.");

  // Check block consistency
  std::map<int, CellType> block_id_map;
  for (const auto& cell : grid->local_cells)
  {
    const int mat_id = cell.material_id;
    if (block_id_map.count(mat_id) == 0)
      block_id_map[mat_id] = cell.GetSubType();
    else
    {
      if (cell.GetSubType() != block_id_map.at(mat_id))
        throw std::logic_error(fname + ": Material id " + std::to_string(mat_id) +
                               " appearing for more than one cell type.");
    }
  }

  // Create unstructured meshes for each material-type pair
  vtkNew<vtkMultiBlockDataSet> grid_blocks;
  int max_dimension = 0;
  {
    vtkNew<vtkUnstructuredGrid> ugrid;
    vtkNew<vtkPoints> points;

    points->SetDataType(VTK_DOUBLE);

    vtkNew<vtkIdTypeArray> global_node_id_list;
    global_node_id_list->SetName("GlobalNodeId");

    vtkNew<vtkIdTypeArray> global_elem_id_list;
    global_elem_id_list->SetName("GlobalElementId");

    vtkNew<vtkIntArray> block_id_list;
    block_id_list->SetName("BlockID");

    // Load vertices
    std::vector<uint64_t> vertex_map(grid->GetGlobalVertexCount(), 0);
    const size_t num_verts = grid->GetGlobalVertexCount();
    for (size_t v = 0; v < num_verts; ++v)
    {
      vertex_map[v] = v;
      const auto& vertex = grid->vertices[v];
      points->InsertNextPoint(vertex.x, vertex.y, vertex.z);

      // Exodus node- and cell indices are 1-based therefore we add a 1 here.
      global_node_id_list->InsertNextValue(static_cast<vtkIdType>(v + 1));
    }

    // Load cells
    for (const auto& cell : grid->local_cells)
    {
      if (cell.GetSubType() == CellType::POLYGON or cell.GetSubType() == CellType::POLYHEDRON)
        throw std::logic_error(fname + ": Cell-subtype \"" + CellTypeName(cell.GetSubType()) +
                               "\" encountered that is not supported by ExodusII.");
      UploadCellGeometryContinuous(cell, vertex_map, ugrid);
      block_id_list->InsertNextValue(cell.material_id);
      max_dimension = std::max(max_dimension, MeshContinuum::GetCellDimension(cell));

      // Exodus node- and cell indices are 1-based therefore we add a 1 here.
      global_elem_id_list->InsertNextValue(static_cast<vtkIdType>(cell.global_id + 1));
    } // for local cells

    // Set arrays
    ugrid->SetPoints(points);
    ugrid->GetPointData()->AddArray(global_node_id_list);
    ugrid->GetCellData()->AddArray(global_elem_id_list);
    ugrid->GetCellData()->AddArray(block_id_list);

    ugrid->GetPointData()->SetActiveGlobalIds("GlobalNodeId");
    ugrid->GetCellData()->SetActiveGlobalIds("GlobalElementId");

    // Set block
    grid_blocks->SetBlock(0, ugrid);

    log.Log() << "Writing grid block "
              << " Number of cells: " << ugrid->GetNumberOfCells()
              << " Number of points: " << ugrid->GetNumberOfPoints();
  } // end of grid_blocks assignment

  // Separate faces by boundary id
  struct FaceInfo
  {
    CellFace const* face_ptr;
    uint64_t source_cell_id;
    int source_face_id;
  };
  std::map<uint64_t, std::vector<FaceInfo>> boundary_id_faces_map;
  for (const auto& cell : grid->local_cells)
  {
    // Here we build a face mapping because OpenSn's face orientation for prisms (wedges) and
    // hexahedrons differ from that of VTK. OpenSn's orientation for prisms and hexes actually
    // matches that of Exodus but VTK assumes the incoming mesh to be conformant to VTK and
    // therefore, internally performs a mapping. Fortunately, the only relevant cell-types, for
    // which a special mapping is required, are the prisms and hexes.
    const size_t num_faces = cell.faces.size();
    std::vector<int> face_mapping(num_faces, 0);
    if (cell.GetSubType() == CellType::WEDGE)
      face_mapping = {2, 3, 4, 0, 1};
    else if (cell.GetSubType() == CellType::HEXAHEDRON)
      face_mapping = {2, 1, 3, 0, 4, 5};
    else
    {
      for (size_t f = 0; f < cell.faces.size(); ++f)
        face_mapping[f] = static_cast<int>(f);
    }

    // Here we store face information as a triplet, i.e., a face pointer, the id of the cell owning
    // it, and the local face index (relative to the cell) of the face.
    int f = 0;
    for (const auto& face : cell.faces)
    {
      if (not face.has_neighbor)
        boundary_id_faces_map[face.neighbor_id].push_back({&face, cell.global_id, face_mapping[f]});
      ++f;
    }
  }

  // Make NodeSets and/or SideSets
  vtkNew<vtkMultiBlockDataSet> nodesets_blocks;
  vtkNew<vtkMultiBlockDataSet> sidesets_blocks;
  for (const auto& [bndry_id, face_list] : boundary_id_faces_map)
  {
    const std::string block_name = grid->GetBoundaryIDMap().at(bndry_id);
    log.Log0Verbose1() << "bid: " + std::to_string(bndry_id) + " name=\"" + block_name + "\"";

    // NodeSet
    {
      vtkNew<vtkUnstructuredGrid> ugrid;
      vtkNew<vtkPoints> points;

      points->SetDataType(VTK_DOUBLE);

      vtkNew<vtkIdTypeArray> node_global_ids;
      node_global_ids->SetName("GlobalNodeId");

      // Build vertex set
      std::set<uint64_t> vid_set;
      for (const auto& face_info : face_list)
        for (uint64_t vid : face_info.face_ptr->vertex_ids)
          vid_set.insert(vid);

      // Build vertex map
      std::vector<uint64_t> vertex_map(grid->GetGlobalVertexCount(), 0);
      {
        uint64_t mapped_id = 0;
        for (uint64_t vid : vid_set)
          vertex_map[vid] = mapped_id++;
      }

      // Load vertices
      for (uint64_t vid : vid_set)
      {
        const auto& vertex = grid->vertices[vid];
        points->InsertNextPoint(vertex.x, vertex.y, vertex.z);

        // Exodus node- and cell indices are 1-based therefore we add a 1 here.
        node_global_ids->InsertNextValue(static_cast<vtkIdType>(vid + 1));
      }

      // Load cells
      for (uint64_t vid : vid_set)
      {
        std::vector<vtkIdType> cell_vids = {static_cast<vtkIdType>(vertex_map[vid])};
        ugrid->InsertNextCell(VTK_VERTEX, static_cast<vtkIdType>(1), cell_vids.data());
      }

      ugrid->SetPoints(points);
      ugrid->GetPointData()->AddArray(node_global_ids);

      ugrid->GetPointData()->SetActiveGlobalIds("GlobalNodeId");

      nodesets_blocks->SetBlock(bndry_id, ugrid);
      nodesets_blocks->GetMetaData(bndry_id)->Set(vtkCompositeDataSet::NAME(), block_name);

      log.Log() << "Writing node set block " << block_name
                << " Number of cells: " << ugrid->GetNumberOfCells()
                << " Number of points: " << ugrid->GetNumberOfPoints();
    }

    // SideSet
    {
      vtkNew<vtkUnstructuredGrid> ugrid;
      vtkNew<vtkPoints> points;

      points->SetDataType(VTK_DOUBLE);

      vtkNew<vtkIdTypeArray> src_cell_global_ids;
      vtkNew<vtkIntArray> src_cell_face_id;
      src_cell_global_ids->SetName("SourceElementId");
      src_cell_face_id->SetName("SourceElementSide");

      // Build vertex set
      std::set<uint64_t> vid_set;
      for (const auto& face_info : face_list)
        for (uint64_t vid : face_info.face_ptr->vertex_ids)
          vid_set.insert(vid);

      // Build vertex map
      std::vector<uint64_t> vertex_map(grid->GetGlobalVertexCount(), 0);
      {
        uint64_t mapped_id = 0;
        for (uint64_t vid : vid_set)
          vertex_map[vid] = mapped_id++;
      }

      // Load vertices
      for (uint64_t vid : vid_set)
      {
        const auto& vertex = grid->vertices[vid];
        points->InsertNextPoint(vertex.x, vertex.y, vertex.z);
      }

      // Load faces
      for (const auto& face_info : face_list)
      {
        UploadFaceGeometry(*face_info.face_ptr, vertex_map, ugrid);
        src_cell_global_ids->InsertNextValue(static_cast<vtkIdType>(face_info.source_cell_id));
        src_cell_face_id->InsertNextValue(face_info.source_face_id);
      }

      ugrid->SetPoints(points);
      ugrid->GetCellData()->AddArray(src_cell_global_ids);
      ugrid->GetCellData()->AddArray(src_cell_face_id);

      sidesets_blocks->SetBlock(bndry_id, ugrid);
      sidesets_blocks->GetMetaData(bndry_id)->Set(vtkCompositeDataSet::NAME(), block_name);

      log.Log() << "Writing side set block " << block_name
                << " Number of cells: " << ugrid->GetNumberOfCells()
                << " Number of points: " << ugrid->GetNumberOfPoints();
    } // End of side-set
  }

  // Write the file
  unsigned int next_block = 0;
  vtkNew<vtkMultiBlockDataSet> main_block;
  main_block->SetBlock(next_block++, grid_blocks);
  if (write_node_sets)
  {
    log.Log0Verbose1() << "Exporting node set";
    main_block->SetBlock(next_block, nodesets_blocks);
    main_block->GetMetaData(next_block++)->Set(vtkCompositeDataSet::NAME(), "Node Sets");
  }
  if (write_side_sets)
  {
    log.Log0Verbose1() << "Exporting side set";
    main_block->SetBlock(next_block, sidesets_blocks);
    main_block->GetMetaData(next_block++)->Set(vtkCompositeDataSet::NAME(), "Side Sets");
  }

  vtkNew<vtkExodusIIWriter> writer;
  writer->SetBlockIdArrayName("BlockID");

  writer->SetFileName(file_name.c_str());
  writer->SetStoreDoubles(1);

  writer->SetInputData(main_block);

  writer->WriteOutGlobalNodeIdArrayOff();
  writer->WriteOutGlobalElementIdArrayOff();
  writer->WriteOutBlockIdArrayOff();

  writer->Write();

  auto em = writer->GetModelMetadata();

  log.Log() << "Num blocks   :  " << em->GetNumberOfBlocks();
  log.Log() << "Num node sets:  " << em->GetNumberOfNodeSets();
  log.Log() << "Num side sets:  " << em->GetNumberOfSideSets();
  log.Log() << "Dimension    :  " << em->GetDimension();
  log.Log() << "Done exporting mesh to ExodusII.";

  opensn::mpi_comm.barrier();
}

void
MeshIO::ToPVTU(const std::shared_ptr<MeshContinuum>& grid, const std::string& file_base_name)
{
  log.Log() << "Exporting mesh to VTK files with base " << file_base_name;

  auto ugrid = PrepareVtkUnstructuredGrid(*grid, false);

  WritePVTUFiles(ugrid, file_base_name);

  log.Log() << "Done exporting mesh to VTK.";

  opensn::mpi_comm.barrier();
}

} // namespace opensn
