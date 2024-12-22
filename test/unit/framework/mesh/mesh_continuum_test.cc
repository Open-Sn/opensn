#include "test/unit/opensn_unit_test.h"
#include "framework/mesh/mesh_continuum/mesh_continuum.h"
#include "framework/mesh/mesh_generator/orthogonal_mesh_generator.h"

using namespace opensn;

class MeshContinuumTest : public OpenSnUnitTest
{
};

/// Helper for building a MeshContinuum for an orthogonal mesh
/// given an array of nodes (one for each dimension)
std::shared_ptr<MeshContinuum>
BuildOrthogonalMesh(const std::vector<std::vector<double>>& node_sets)
{
  ParameterBlock array("node_sets");
  for (std::size_t i = 0; i < node_sets.size(); ++i)
    array.AddParameter(ParameterBlock(std::to_string(i + 1), node_sets[i]));
  array.ChangeToArray();

  ParameterBlock block("");
  block.AddParameter(array);

  auto params = OrthogonalMeshGenerator::GetInputParameters();
  params.AssignParameters(block);
  auto grid_ptr = OrthogonalMeshGenerator(params).Execute();
  return grid_ptr;
}

/// Helper for the PointInsideCellXD tests
void
TestPointInsideCell(const MeshContinuum& grid)
{
  // Centroid is contained within cell whose centroid it is
  for (const auto& cell : grid.local_cells)
    for (const auto& other_cell : grid.local_cells)
    {
      const auto same_cell = cell.global_id == other_cell.global_id;
      const auto within = grid.CheckPointInsideCell(other_cell, cell.centroid);
      EXPECT_EQ(same_cell, within);
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
        EXPECT_EQ(has_vertex, within);
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
        EXPECT_EQ(same_cell_or_neighbor, within);
      }

  if (grid.GetDimension() > 1)
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

  if (grid.GetDimension() == 3)
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
            EXPECT_EQ(same_cell, within);
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
              EXPECT_EQ(same_cell_or_neighbor, within);
            }
          }
        }
  }
}

TEST_F(MeshContinuumTest, PointInsideCell1D)
{
  const auto grid_ptr = BuildOrthogonalMesh({{-1.0, -0.75, 0.0, 1.0, 2.0}});
  TestPointInsideCell(*grid_ptr);
}

TEST_F(MeshContinuumTest, PointInsideCell2D)
{
  const auto grid_ptr = BuildOrthogonalMesh({{-1.0, -0.75, 0.0, 1.0}, {0.0, 0.5, 1.0}});
  TestPointInsideCell(*grid_ptr);
}

TEST_F(MeshContinuumTest, PointInsideCell3D)
{
  const auto grid_ptr = BuildOrthogonalMesh({{-1.0, 1.0}, {0.0, 0.5, 1.0}, {-1.0, 0.0, 1.0}});
  TestPointInsideCell(*grid_ptr);
}

/// Helper for the PointInsideCellFaceXD tests
void
TestPointInsideCellFace(const MeshContinuum& grid)
{
  // Vertices contained within faces that have those vertices
  for (uint64_t vi = 0; vi < grid.GetGlobalVertexCount(); ++vi)
    for (const auto& cell : grid.local_cells)
      for (std::size_t face_i = 0; face_i < cell.faces.size(); ++face_i)
      {
        const auto& face = cell.faces[face_i];
        const auto has_vertex =
          std::find(face.vertex_ids.begin(), face.vertex_ids.end(), vi) != face.vertex_ids.end();
        const auto within = grid.CheckPointInsideCellFace(cell, face_i, grid.vertices[vi]);
        EXPECT_EQ(has_vertex, within);
      }

  // Cell centroids not contained within any faces
  for (const auto& cell : grid.local_cells)
    for (const auto& other_cell : grid.local_cells)
      for (std::size_t other_face_i = 0; other_face_i < other_cell.faces.size(); ++other_face_i)
        EXPECT_FALSE(grid.CheckPointInsideCellFace(other_cell, other_face_i, cell.centroid));

  if (grid.GetDimension() > 1)
  {
    // Face centroids contained within faces
    for (const auto& cell : grid.local_cells)
      for (const auto& face : cell.faces)
        for (const auto& other_cell : grid.local_cells)
          for (std::size_t other_face_i = 0; other_face_i < other_cell.faces.size(); ++other_face_i)
          {
            const auto& other_face = other_cell.faces[other_face_i];
            const auto same_face = other_face.centroid.AbsoluteEquals(face.centroid);
            const auto within =
              grid.CheckPointInsideCellFace(other_cell, other_face_i, face.centroid);
            EXPECT_EQ(same_face, within);
          }
  }

  if (grid.GetDimension() > 2)
  {
    // Face edge centroids contained within faces that have the vertices that contain the edge
    for (const auto& cell : grid.local_cells)
      for (const auto& face : cell.faces)
        for (std::size_t i = 0; i < face.vertex_ids.size(); ++i)
        {
          const auto vi1 = face.vertex_ids[i];
          const auto vi2 =
            (i == face.vertex_ids.size() - 1) ? face.vertex_ids[0] : face.vertex_ids[i + 1];
          const auto edge_centroid = (grid.vertices[vi1] + grid.vertices[vi2]) / 2;
          for (const auto& other_cell : grid.local_cells)
            for (std::size_t other_face_i = 0; other_face_i < other_cell.faces.size();
                 ++other_face_i)
            {
              const auto& other_face = other_cell.faces[other_face_i];
              const auto has_edges =
                (std::find(other_face.vertex_ids.begin(), other_face.vertex_ids.end(), vi1) !=
                 other_face.vertex_ids.end()) &&
                (std::find(other_face.vertex_ids.begin(), other_face.vertex_ids.end(), vi2) !=
                 other_face.vertex_ids.end());
              const auto within =
                grid.CheckPointInsideCellFace(other_cell, other_face_i, edge_centroid);
              EXPECT_EQ(has_edges, within);
            }
        }
  }
}

TEST_F(MeshContinuumTest, PointInsideCellFace1D)
{
  const auto grid_ptr = BuildOrthogonalMesh({{-1.0, 0.0, 1.0, 2.0}});
  TestPointInsideCellFace(*grid_ptr);
}

TEST_F(MeshContinuumTest, PointInsideCellFace2D)
{
  const auto grid_ptr = BuildOrthogonalMesh({{-1.0, -0.75, 0.0}, {0.0, 0.5, 1.0}});
  TestPointInsideCellFace(*grid_ptr);
}

TEST_F(MeshContinuumTest, PointInsideCellFace3D)
{
  const auto grid_ptr = BuildOrthogonalMesh({{-1.0, 1.0}, {0.0, 0.5, 1.0}, {-1.0, 0.0, 1.0}});
  TestPointInsideCellFace(*grid_ptr);
}
