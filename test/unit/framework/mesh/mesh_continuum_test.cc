#include "gtest/gtest.h"
#include "test/unit/common/mesh_builders.h"
#include "framework/mesh/mesh/mesh.h"
#include "framework/runtime.h"

using namespace opensn;

/// Helper for the PointInsideCellXD tests
void
TestPointInsideCell(const std::shared_ptr<Mesh> grid)
{
  // Centroid is contained within cell whose centroid it is
  for (std::size_t cell_local_id = 0; cell_local_id < grid->GetLocalCellCount(); ++cell_local_id)
  {
    const auto& cell = grid->GetLocalCell(cell_local_id);
    for (std::size_t other_cell_local_id = 0; other_cell_local_id < grid->GetLocalCellCount();
         ++other_cell_local_id)
    {
      const auto& other_cell = grid->GetLocalCell(other_cell_local_id);
      const auto same_cell = cell.global_id == other_cell.global_id;
      const auto within = grid->CheckPointInsideCell(other_cell_local_id, cell.centroid);
      EXPECT_EQ(same_cell, within);
    }
  }

  // Vertices are contained within cells
  for (std::size_t cell_local_id = 0; cell_local_id < grid->GetLocalCellCount(); ++cell_local_id)
  {
    const auto& cell_vertex_ids = grid->GetCellConnectivity(cell_local_id);
    for (const auto vi : cell_vertex_ids)
      for (std::size_t other_cell_local_id = 0; other_cell_local_id < grid->GetLocalCellCount();
           ++other_cell_local_id)
      {
        const auto& other_cell_vertex_ids = grid->GetCellConnectivity(other_cell_local_id);
        const auto has_vertex =
          std::find(other_cell_vertex_ids.begin(), other_cell_vertex_ids.end(), vi) !=
          other_cell_vertex_ids.end();
        const auto within = grid->CheckPointInsideCell(other_cell_local_id, grid->GlobalVertex(vi));
        EXPECT_EQ(has_vertex, within);
      }
  }

  // Face centroids are contained within cells (including neighbors)
  for (const auto& cell : grid->GetLocalCells())
    for (const auto& face : cell.faces)
      for (std::size_t other_cell_local_id = 0; other_cell_local_id < grid->GetLocalCellCount();
           ++other_cell_local_id)
      {
        const auto& other_cell = grid->GetLocalCell(other_cell_local_id);
        const auto same_cell_or_neighbor =
          (cell.global_id == other_cell.global_id ||
           (face.has_neighbor && face.neighbor_id == other_cell.global_id));
        const auto within = grid->CheckPointInsideCell(other_cell_local_id, face.centroid);
        EXPECT_EQ(same_cell_or_neighbor, within);
      }

  if (grid->GetDimension() > 1)
  {
    // Face edge centroids contained
    for (std::size_t cell_local_id = 0; cell_local_id < grid->GetLocalCellCount(); ++cell_local_id)
    {
      const auto& cell = grid->GetLocalCell(cell_local_id);
      for (size_t f = 0; f < cell.faces.size(); ++f)
      {
        const auto& face = cell.faces[f];
        const auto face_vertex_ids = grid->GetCellFaceConnectivity(cell_local_id, f);
        for (size_t side = 0; side < face_vertex_ids.size(); ++side)
        {
          const size_t sp1 = (side < (face_vertex_ids.size() - 1)) ? side + 1 : 0;
          const auto& v0 = grid->GlobalVertex(face_vertex_ids[side]);
          const auto& v1 = grid->GlobalVertex(face_vertex_ids[sp1]);
          const auto c = (v0 + v1) / 2.0;
          EXPECT_TRUE(grid->CheckPointInsideCell(cell_local_id, c));
        }
      }
    }
  }

  if (grid->GetDimension() == 3)
  {
    // Tetrahedral centroids contained
    for (const auto& cell : grid->GetLocalCells())
    {
      const auto cell_local_id = grid->MapCellGlobalID2LocalID(cell.global_id);
      for (size_t f = 0; f < cell.faces.size(); ++f)
      {
        const auto& face = cell.faces[f];
        const auto face_vertex_ids = grid->GetCellFaceConnectivity(cell_local_id, f);
        for (size_t side = 0; side < face_vertex_ids.size(); ++side)
        {
          const size_t sp1 = (side < (face_vertex_ids.size() - 1)) ? side + 1 : 0;
          const auto& v0 = grid->GlobalVertex(face_vertex_ids[side]);
          const auto& v1 = face.centroid;
          const auto& v2 = grid->GlobalVertex(face_vertex_ids[sp1]);
          const auto& v3 = cell.centroid;
          const auto c = (v0 + v1 + v2 + v3) / 4.0;
          for (std::size_t other_cell_local_id = 0; other_cell_local_id < grid->GetLocalCellCount();
               ++other_cell_local_id)
          {
            const auto& other_cell = grid->GetLocalCell(other_cell_local_id);
            const auto same_cell = cell.global_id == other_cell.global_id;
            const auto within = grid->CheckPointInsideCell(other_cell_local_id, c);
            EXPECT_EQ(same_cell, within);
          }
        }
      }
    }

    // Tetrahedral face centroids contained
    for (const auto& cell : grid->GetLocalCells())
    {
      const auto cell_local_id = grid->MapCellGlobalID2LocalID(cell.global_id);
      for (size_t f = 0; f < cell.faces.size(); ++f)
      {
        const auto& face = cell.faces[f];
        const auto face_vertex_ids = grid->GetCellFaceConnectivity(cell_local_id, f);
        for (size_t side = 0; side < face_vertex_ids.size(); ++side)
        {
          const auto tet_face_vertices = grid->GetTetrahedralFaceVertices(cell_local_id, f, side);
          for (const auto& v : tet_face_vertices)
            for (std::size_t other_cell_local_id = 0;
                 other_cell_local_id < grid->GetLocalCellCount();
                 ++other_cell_local_id)
            {
              const auto& other_cell = grid->GetLocalCell(other_cell_local_id);
              const auto same_cell_or_neighbor =
                (cell.global_id == other_cell.global_id ||
                 (face.has_neighbor && face.neighbor_id == other_cell.global_id));
              const auto within = grid->CheckPointInsideCell(other_cell_local_id, face.centroid);
              EXPECT_EQ(same_cell_or_neighbor, within);
            }
        }
      }
    }
  }
}

TEST(MeshTest, PointInsideCell1D)
{
  const auto grid_ptr = BuildOrthogonalMesh({{-1.0, -0.75, 0.0, 1.0, 2.0}});
  TestPointInsideCell(grid_ptr);
}

TEST(MeshTest, PointInsideCell2D)
{
  const auto grid_ptr = BuildOrthogonalMesh({{-1.0, -0.75, 0.0, 1.0}, {0.0, 0.5, 1.0}});
  TestPointInsideCell(grid_ptr);
}

TEST(MeshTest, PointInsideCell3D)
{
  const auto grid_ptr = BuildOrthogonalMesh({{-1.0, 1.0}, {0.0, 0.5, 1.0}, {-1.0, 0.0, 1.0}});
  TestPointInsideCell(grid_ptr);
}

/// Helper for the PointInsideCellFaceXD tests
void
TestPointInsideCellFace(const std::shared_ptr<Mesh> grid)
{
  // Vertices contained within faces that have those vertices
  for (uint64_t vi = 0; vi < grid->GetGlobalVertexCount(); ++vi)
    for (const auto& cell : grid->GetLocalCells())
    {
      const auto cell_local_id = grid->MapCellGlobalID2LocalID(cell.global_id);
      for (std::size_t face_i = 0; face_i < cell.faces.size(); ++face_i)
      {
        const auto face_vertex_ids = grid->GetCellFaceConnectivity(cell_local_id, face_i);
        const auto has_vertex =
          std::find(face_vertex_ids.begin(), face_vertex_ids.end(), vi) != face_vertex_ids.end();
        const auto within = grid->CheckPointInsideCellFace(cell_local_id, face_i, grid->GlobalVertex(vi));
        EXPECT_EQ(has_vertex, within);
      }
    }

  // Cell centroids not contained within any faces
  for (const auto& cell : grid->GetLocalCells())
    for (const auto& other_cell : grid->GetLocalCells())
    {
      const auto other_cell_local_id = grid->MapCellGlobalID2LocalID(other_cell.global_id);
      for (std::size_t other_face_i = 0; other_face_i < other_cell.faces.size(); ++other_face_i)
        EXPECT_FALSE(grid->CheckPointInsideCellFace(other_cell_local_id, other_face_i, cell.centroid));
    }

  if (grid->GetDimension() > 1)
  {
    // Face centroids contained within faces
    for (const auto& cell : grid->GetLocalCells())
      for (const auto& face : cell.faces)
        for (const auto& other_cell : grid->GetLocalCells())
        {
          const auto other_cell_local_id = grid->MapCellGlobalID2LocalID(other_cell.global_id);
          for (std::size_t other_face_i = 0; other_face_i < other_cell.faces.size(); ++other_face_i)
          {
            const auto& other_face = other_cell.faces[other_face_i];
            const auto same_face = other_face.centroid.AbsoluteEquals(face.centroid);
            const auto within =
              grid->CheckPointInsideCellFace(other_cell_local_id, other_face_i, face.centroid);
            EXPECT_EQ(same_face, within);
          }
        }
  }

  if (grid->GetDimension() > 2)
  {
    // Face edge centroids contained within faces that have the vertices that contain the edge
    for (const auto& cell : grid->GetLocalCells())
    {
      const auto cell_local_id = grid->MapCellGlobalID2LocalID(cell.global_id);
      for (std::size_t face_i = 0; face_i < cell.faces.size(); ++face_i)
      {
        const auto face_vertex_ids = grid->GetCellFaceConnectivity(cell_local_id, face_i);
        for (std::size_t i = 0; i < face_vertex_ids.size(); ++i)
        {
          const auto vi1 = face_vertex_ids[i];
          const auto vi2 =
            (i == face_vertex_ids.size() - 1) ? face_vertex_ids[0] : face_vertex_ids[i + 1];
          const auto edge_centroid = (grid->GlobalVertex(vi1) + grid->GlobalVertex(vi2)) / 2;
          for (const auto& other_cell : grid->GetLocalCells())
          {
            const auto other_cell_local_id = grid->MapCellGlobalID2LocalID(other_cell.global_id);
            for (std::size_t other_face_i = 0; other_face_i < other_cell.faces.size();
                 ++other_face_i)
            {
              const auto& other_face = other_cell.faces[other_face_i];
              const auto other_face_vertex_ids = grid->GetCellFaceConnectivity(other_cell_local_id, other_face_i);
              const auto has_edges =
                (std::find(other_face_vertex_ids.begin(), other_face_vertex_ids.end(), vi1) !=
                 other_face_vertex_ids.end()) &&
                (std::find(other_face_vertex_ids.begin(), other_face_vertex_ids.end(), vi2) !=
                 other_face_vertex_ids.end());
              const auto within =
                grid->CheckPointInsideCellFace(other_cell_local_id, other_face_i, edge_centroid);
              EXPECT_EQ(has_edges, within);
            }
          }
        }
      }
    }
  }
}

TEST(MeshTest, PointInsideCellFace1D)
{
  if (opensn::mpi_comm.size() != 1)
    return;

  const auto grid_ptr = BuildOrthogonalMesh({{-1.0, 0.0, 1.0, 2.0}});
  TestPointInsideCellFace(grid_ptr);
}

TEST(MeshTest, PointInsideCellFace2D)
{
  const auto grid_ptr = BuildOrthogonalMesh({{-1.0, -0.75, 0.0}, {0.0, 0.5, 1.0}});
  TestPointInsideCellFace(grid_ptr);
}

TEST(MeshTest, PointInsideCellFace3D)
{
  const auto grid_ptr = BuildOrthogonalMesh({{-1.0, 1.0}, {0.0, 0.5, 1.0}, {-1.0, 0.0, 1.0}});
  TestPointInsideCellFace(grid_ptr);
}
