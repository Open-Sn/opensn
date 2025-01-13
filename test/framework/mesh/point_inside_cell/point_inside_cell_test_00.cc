#include "test/src/unit_tester.h"

#include "framework/mesh/mesh.h"
#include "framework/mesh/mesh_continuum/mesh_continuum.h"

using namespace opensn;

namespace unit_tests
{

ParameterBlock point_inside_cell_test_00(const InputParameters& params);

RegisterWrapperFunctionInNamespace(unit_tests,
                                   point_inside_cell_test_00,
                                   nullptr,
                                   point_inside_cell_test_00);

/// Exhaustive testing for MeshContinuum::CheckPointInsideCell
ParameterBlock
point_inside_cell_test_00(const InputParameters&)
{
  const auto grid_ptr = GetCurrentMesh();
  const auto& grid = *grid_ptr;

  SETUP_UNIT_TESTER();

  // Centroid is contained within cell whose centroid it is
  for (const auto& cell : grid.local_cells)
    for (const auto& other_cell : grid.local_cells)
    {
      const auto same_cell = cell.global_id == other_cell.global_id;
      const auto within = grid.CheckPointInsideCell(other_cell, cell.centroid);
      EXPECT_EQUAL(same_cell, within);
    }

  // Vertices are contained within cells
  for (const auto& cell : grid.local_cells)
    for (const auto vi : cell.vertex_ids)
      for (const auto& other_cell : grid.local_cells)
      {
        const auto has_vertex =
          std::find(other_cell.vertex_ids.begin(), other_cell.vertex_ids.end(), vi) !=
          other_cell.vertex_ids.end();
        const auto within = grid.CheckPointInsideCell(other_cell, grid.vertices[vi]);
        EXPECT_EQUAL(has_vertex, within);
      }

  // Face centroids are contained within cells (including neighbors)
  for (const auto& cell : grid.local_cells)
    for (const auto& face : cell.faces)
      for (const auto& other_cell : grid.local_cells)
      {
        const auto same_cell_or_neighbor =
          (cell.global_id == other_cell.global_id ||
           (face.has_neighbor && face.neighbor_id == other_cell.global_id));
        const auto within = grid.CheckPointInsideCell(other_cell, face.centroid);
        EXPECT_EQUAL(same_cell_or_neighbor, within);
      }

  if (grid.Dimension() > 1)
  {
    // Face edge centroids contained
    for (const auto& cell : grid.local_cells)
      for (const auto& face : cell.faces)
        for (size_t side = 0; side < face.vertex_ids.size(); ++side)
        {
          const size_t sp1 = (side < (face.vertex_ids.size() - 1)) ? side + 1 : 0;
          const auto& v0 = grid.vertices[face.vertex_ids[side]];
          const auto& v1 = grid.vertices[face.vertex_ids[sp1]];
          const auto c = (v0 + v1) / 2.0;
          EXPECT_TRUE(grid.CheckPointInsideCell(cell, c));
        }
  }

  if (grid.Dimension() == 3)
  {
    // Tetrahedral centroids contained
    for (const auto& cell : grid.local_cells)
      for (const auto& face : cell.faces)
        for (size_t side = 0; side < face.vertex_ids.size(); ++side)
        {
          const size_t sp1 = (side < (face.vertex_ids.size() - 1)) ? side + 1 : 0;
          const auto& v0 = grid.vertices[face.vertex_ids[side]];
          const auto& v1 = face.centroid;
          const auto& v2 = grid.vertices[face.vertex_ids[sp1]];
          const auto& v3 = cell.centroid;
          const auto c = (v0 + v1 + v2 + v3) / 4.0;
          for (const auto& other_cell : grid.local_cells)
          {
            const auto same_cell = cell.global_id == other_cell.global_id;
            const auto within = grid.CheckPointInsideCell(other_cell, c);
            EXPECT_EQUAL(same_cell, within);
          }
        }

    // Tetrahedral face centroids contained
    for (const auto& cell : grid.local_cells)
      for (const auto& face : cell.faces)
        for (size_t side = 0; side < face.vertex_ids.size(); ++side)
        {
          const auto tet_face_vertices = grid.GetTetrahedralFaceVertices(cell, face, side);
          for (const auto& v : tet_face_vertices)
          {
            const auto c = (v[0] + v[1] + v[2]) / 3.0;
            for (const auto& other_cell : grid.local_cells)
            {
              const auto same_cell_or_neighbor =
                (cell.global_id == other_cell.global_id ||
                 (face.has_neighbor && face.neighbor_id == other_cell.global_id));
              const auto within = grid.CheckPointInsideCell(other_cell, face.centroid);
              EXPECT_EQUAL(same_cell_or_neighbor, within);
            }
          }
        }
  }

  FINALIZE_UNIT_TESTER();

  return ParameterBlock();
}

} //  namespace unit_tests
