#include "chi_unpartitioned_mesh.h"

#include <vtkPolygon.h>
#include <vtkLine.h>
#include <vtkVertex.h>

// ###################################################################
/**Creates a raw polyhedron cell from a vtk-polyhedron.*/
chi_mesh::UnpartitionedMesh::LightWeightCell*
chi_mesh::UnpartitionedMesh::CreateCellFromVTKPolyhedron(vtkCell* vtk_cell)
{
  const std::string fname =
    "chi_mesh::UnpartitionedMesh::CreateCellFromVTKPolyhedron";

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
      std::vector<std::vector<uint64_t>> face_vids = {{1, 2, 6, 5},
                                                      {3, 0, 4, 7},
                                                      {2, 3, 7, 6},
                                                      {0, 1, 5, 4},
                                                      {4, 5, 6, 7},
                                                      {3, 2, 1, 0}};
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
      std::vector<std::vector<uint64_t>> face_vids = {{0, 2, 1},
                                                      {0, 1, 3},
                                                      {0, 3, 2},
                                                      {3, 1, 2}};
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

// ###################################################################
/**Creates a raw polygon cell from a vtk-polygon.*/
chi_mesh::UnpartitionedMesh::LightWeightCell*
chi_mesh::UnpartitionedMesh::CreateCellFromVTKPolygon(vtkCell* vtk_cell)
{
  const std::string fname =
    "chi_mesh::UnpartitionedMesh::CreateCellFromVTKPolygon";

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
    auto v1_id = (f < (num_cfaces - 1)) ? poly_cell->vertex_ids[f + 1]
                                        : poly_cell->vertex_ids[0];

    face.vertex_ids.reserve(2);
    face.vertex_ids.push_back(v0_id);
    face.vertex_ids.push_back(v1_id);

    poly_cell->faces.push_back(face);
  }

  return poly_cell;
}

// ###################################################################
/**Creates a raw slab cell from a vtk-line.*/
chi_mesh::UnpartitionedMesh::LightWeightCell*
chi_mesh::UnpartitionedMesh::CreateCellFromVTKLine(vtkCell* vtk_cell)
{
  const std::string fname =
    "chi_mesh::UnpartitionedMesh::CreateCellFromVTKPolygon";

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

// ###################################################################
/**Creates a raw point cell from a vtk-vertex.*/
chi_mesh::UnpartitionedMesh::LightWeightCell*
chi_mesh::UnpartitionedMesh::CreateCellFromVTKVertex(vtkCell* vtk_cell)
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