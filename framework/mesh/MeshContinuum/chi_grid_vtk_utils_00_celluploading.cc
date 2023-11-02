#include "chi_grid_vtk_utils.h"

#include "chi_meshcontinuum.h"

#include <vtkPoints.h>
#include <vtkUnstructuredGrid.h>

//###################################################################
/**Uploads vertices and cells to an unstructured grid.*/
void
chi_mesh::UploadCellGeometryDiscontinuous(const chi_mesh::MeshContinuum& grid,
                                          const chi_mesh::Cell& cell,
                                          int64_t& node_counter,
                                          vtkNew<vtkPoints>& points,
                                          vtkNew<vtkUnstructuredGrid>& ugrid)
{
  size_t num_verts = cell.vertex_ids_.size();

  std::vector<vtkIdType> cell_vids(num_verts);
  for (size_t v = 0; v < num_verts; v++)
  {
    uint64_t vgi = cell.vertex_ids_[v];
    std::vector<double> d_node(3);
    d_node[0] = grid.vertices[vgi].x;
    d_node[1] = grid.vertices[vgi].y;
    d_node[2] = grid.vertices[vgi].z;

    points->InsertPoint(node_counter, d_node.data());
    cell_vids[v] = node_counter++;
  }

  if (cell.Type() == chi_mesh::CellType::SLAB)
  {
    ugrid->InsertNextCell(VTK_LINE, static_cast<vtkIdType>(num_verts), cell_vids.data());
  }
  if (cell.Type() == chi_mesh::CellType::POLYGON)
  {
    int vtk_subtype;
    switch (cell.SubType())
    {
      case CellType::POLYGON:
        vtk_subtype = VTK_POLYGON;
        break;
      case CellType::QUADRILATERAL:
        vtk_subtype = VTK_QUAD;
        break;
      case CellType::TRIANGLE:
        vtk_subtype = VTK_TRIANGLE;
        break;
      default:
        vtk_subtype = VTK_POLYGON;
        break;
    }

    ugrid->InsertNextCell(vtk_subtype, static_cast<vtkIdType>(num_verts), cell_vids.data());
  }
  if (cell.Type() == chi_mesh::CellType::POLYHEDRON)
  {
    // Build polyhedron faces
    std::vector<vtkIdType> faces_vids;

    size_t num_faces = cell.faces_.size();
    for (auto& face : cell.faces_)
    {
      size_t num_fverts = face.vertex_ids_.size();
      std::vector<vtkIdType> face_info(num_fverts);
      for (size_t fv = 0; fv < num_fverts; fv++)
      {
        size_t v = 0;
        for (size_t cv = 0; cv < num_verts; ++cv)
          if (cell.vertex_ids_[cv] == face.vertex_ids_[fv])
          {
            v = cv;
            break;
          }

        face_info[fv] = cell_vids[v];
      }

      faces_vids.push_back(static_cast<vtkIdType>(num_fverts));
      for (auto vid : face_info)
        faces_vids.push_back(vid);
    } // for f

    int vtk_subtype;
    switch (cell.SubType())
    {
      case CellType::POLYHEDRON:
        vtk_subtype = VTK_POLYHEDRON;
        break;
      case CellType::PYRAMID:
        vtk_subtype = VTK_PYRAMID;
        break;
      case CellType::WEDGE:
        vtk_subtype = VTK_WEDGE;
        break;
      case CellType::HEXAHEDRON:
        vtk_subtype = VTK_HEXAHEDRON;
        break;
      case CellType::TETRAHEDRON:
        vtk_subtype = VTK_TETRA;
        break;
      default:
        vtk_subtype = VTK_POLYHEDRON;
        break;
    }

    ugrid->InsertNextCell(vtk_subtype,
                          static_cast<vtkIdType>(num_verts),
                          cell_vids.data(),
                          static_cast<vtkIdType>(num_faces),
                          faces_vids.data());
  } // polyhedron
}

//###################################################################
/**Uploads vertices and cells to an unstructured grid.*/
void
chi_mesh::UploadCellGeometryContinuous(const chi_mesh::Cell& cell,
                                       const std::vector<uint64_t>& vertex_map,
                                       vtkNew<vtkUnstructuredGrid>& ugrid)
{
  size_t num_verts = cell.vertex_ids_.size();

  std::vector<vtkIdType> cell_vids(num_verts);
  for (size_t v = 0; v < num_verts; v++)
    cell_vids[v] = static_cast<vtkIdType>(vertex_map[cell.vertex_ids_[v]]);

  if (cell.Type() == chi_mesh::CellType::SLAB)
  {
    ugrid->InsertNextCell(VTK_LINE, static_cast<vtkIdType>(num_verts), cell_vids.data());
  }
  if (cell.Type() == chi_mesh::CellType::POLYGON)
  {
    int vtk_subtype;
    switch (cell.SubType())
    {
      case CellType::POLYGON:
        vtk_subtype = VTK_POLYGON;
        break;
      case CellType::QUADRILATERAL:
        vtk_subtype = VTK_QUAD;
        break;
      case CellType::TRIANGLE:
        vtk_subtype = VTK_TRIANGLE;
        break;
      default:
        vtk_subtype = VTK_POLYGON;
        break;
    }

    ugrid->InsertNextCell(vtk_subtype, static_cast<vtkIdType>(num_verts), cell_vids.data());
  }
  if (cell.Type() == chi_mesh::CellType::POLYHEDRON)
  {
    // Build polyhedron faces
    std::vector<vtkIdType> faces_vids;

    int vtk_subtype;
    switch (cell.SubType())
    {
      case CellType::POLYHEDRON:
        vtk_subtype = VTK_POLYHEDRON;
        break;
      case CellType::PYRAMID:
        vtk_subtype = VTK_PYRAMID;
        break;
      case CellType::WEDGE:
        vtk_subtype = VTK_WEDGE;
        break;
      case CellType::HEXAHEDRON:
        vtk_subtype = VTK_HEXAHEDRON;
        break;
      case CellType::TETRAHEDRON:
        vtk_subtype = VTK_TETRA;
        break;
      default:
        vtk_subtype = VTK_POLYHEDRON;
        break;
    }

    switch (cell.SubType())
    {
      case CellType::POLYHEDRON:
      {
        size_t num_faces = cell.faces_.size();
        for (auto& face : cell.faces_)
        {
          size_t num_fverts = face.vertex_ids_.size();
          std::vector<vtkIdType> face_info(num_fverts);
          for (size_t fv = 0; fv < num_fverts; fv++)
          {
            size_t v = 0;
            for (size_t cv = 0; cv < num_verts; ++cv)
              if (cell.vertex_ids_[cv] == face.vertex_ids_[fv])
              {
                v = cv;
                break;
              }

            face_info[fv] = cell_vids[v];
          }

          faces_vids.push_back(static_cast<vtkIdType>(num_fverts));
          for (auto vid : face_info)
            faces_vids.push_back(vid);
        } // for f

        ugrid->InsertNextCell(vtk_subtype,
                              static_cast<vtkIdType>(num_verts),
                              cell_vids.data(),
                              static_cast<vtkIdType>(num_faces),
                              faces_vids.data());
        break;
      }
      default:
        ugrid->InsertNextCell(vtk_subtype, static_cast<vtkIdType>(num_verts), cell_vids.data());
    }
  } // polyhedron
}

//###################################################################
/**Uploads vertices and cells to an unstructured grid.*/
void
chi_mesh::UploadFaceGeometry(const chi_mesh::CellFace& cell_face,
                             const std::vector<uint64_t>& vertex_map,
                             vtkNew<vtkUnstructuredGrid>& ugrid)
{
  const size_t num_verts = cell_face.vertex_ids_.size();

  std::vector<vtkIdType> cell_vids;
  for (uint64_t vid : cell_face.vertex_ids_)
    cell_vids.push_back(static_cast<vtkIdType>(vertex_map[vid]));

  if (num_verts == 1)
  {
    ugrid->InsertNextCell(VTK_VERTEX, static_cast<vtkIdType>(num_verts), cell_vids.data());
  }
  if (num_verts == 2)
  {
    ugrid->InsertNextCell(VTK_LINE, static_cast<vtkIdType>(num_verts), cell_vids.data());
  }
  if (num_verts >= 3)
  {
    int vtk_subtype;
    switch (num_verts)
    {
      case 3:
        vtk_subtype = VTK_TRIANGLE;
        break;
      case 4:
        vtk_subtype = VTK_QUAD;
        break;
      default:
        vtk_subtype = VTK_POLYGON;
        break;
    }

    ugrid->InsertNextCell(vtk_subtype, static_cast<vtkIdType>(num_verts), cell_vids.data());
  }
}