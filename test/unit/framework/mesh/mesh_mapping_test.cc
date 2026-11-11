#include "gtest/gtest.h"
#include "test/unit/common/mesh_builders.h"
#include "framework/mesh/mesh/mesh.h"
#include "framework/mesh/mesh_mapping/mesh_mapping.h"
#include "framework/runtime.h"
#include <filesystem>

using namespace opensn;

namespace
{

void
TestMapping(const std::shared_ptr<Mesh> fine_grid, const std::shared_ptr<Mesh> coarse_grid)
{
  MeshMapping mesh_mapping;
  mesh_mapping.Build(fine_grid, coarse_grid);

  // Volumetric mapping; check all coarse cells against all fine cells
  for (std::size_t clid = 0; clid < coarse_grid->GetLocalCellCount(); ++clid)
  {
    const auto& coarse_cell = coarse_grid->GetLocalCell(clid);
    const auto& coarse_mapping = mesh_mapping.GetCoarseMapping(coarse_cell);
    for (std::size_t flid = 0; flid < fine_grid->GetLocalCellCount(); ++flid)
    {
      const auto& fine_cell = fine_grid->GetLocalCell(flid);
      const auto in_coarse_mapping =
        std::find(coarse_mapping.fine_cells.begin(), coarse_mapping.fine_cells.end(), &fine_cell) !=
        coarse_mapping.fine_cells.end();
      const auto& fine_mapping = mesh_mapping.GetFineMapping(fine_cell);
      const auto in_fine_mapping = fine_mapping.coarse_cell == &coarse_cell;
      EXPECT_EQ(in_coarse_mapping, in_fine_mapping);

      const auto fine_within_coarse =
        coarse_grid->CheckPointInsideCell(coarse_cell, fine_cell.centroid);
      EXPECT_EQ(in_coarse_mapping, fine_within_coarse);
    }
  }

  // Surface mapping. Check all fine cell faces against all coarse cell faces
  for (std::size_t flid = 0; flid < fine_grid->GetLocalCellCount(); ++flid)
  {
    const auto& fine_cell = fine_grid->GetLocalCell(flid);
    const auto& fine_mapping = mesh_mapping.GetFineMapping(fine_cell);
    for (std::size_t fine_face_i = 0; fine_face_i < fine_cell.faces.size(); ++fine_face_i)
    {
      const auto& fine_face = fine_cell.faces[fine_face_i];
      bool found_mapping = false;
      for (std::size_t clid = 0; clid < coarse_grid->GetLocalCellCount(); ++clid)
      {
        const auto& coarse_cell = coarse_grid->GetLocalCell(clid);
        const auto& coarse_mapping = mesh_mapping.GetCoarseMapping(coarse_cell);
        for (std::size_t coarse_face_i = 0; coarse_face_i < coarse_cell.faces.size();
             ++coarse_face_i)
        {
          const auto in_fine_mapping = (fine_mapping.coarse_cell == &coarse_cell) &&
                                       (fine_mapping.coarse_faces[fine_face_i] == coarse_face_i);
          const auto& coarse_face_mapping = coarse_mapping.fine_faces[coarse_face_i];
          const auto in_coarse_mapping =
            std::find(coarse_face_mapping.begin(),
                      coarse_face_mapping.end(),
                      std::make_pair(&fine_cell, fine_face_i)) != coarse_face_mapping.end();
          EXPECT_EQ(in_fine_mapping, in_coarse_mapping);

          if (in_fine_mapping)
          {
            EXPECT_FALSE(found_mapping);
            found_mapping = true;
          }
          const auto fine_face_within_coarse =
            coarse_grid->CheckPointInsideCellFace(coarse_cell, coarse_face_i, fine_face.centroid);
          const auto fine_within_coarse =
            coarse_grid->CheckPointInsideCell(coarse_cell, fine_cell.centroid);
          EXPECT_EQ(in_fine_mapping, fine_face_within_coarse && fine_within_coarse);
        }
      }
      if (!found_mapping)
      {
        EXPECT_EQ(MeshMapping::invalid_face_index, fine_mapping.coarse_faces[fine_face_i]);
      }
    }
  }
}

} // namespace

TEST(MeshMappingTest, Test1D)
{
  if (opensn::mpi_comm.size() != 1)
    return;

  const auto fine_grid_ptr = BuildOrthogonalMesh({{-1.0, -0.8, -0.75, -0.5, -0.25, 0.0, 0.5, 1.0}});
  const auto coarse_grid_ptr = BuildOrthogonalMesh({{-1.0, -0.5, 0.0, 0.5, 1.0}});
  TestMapping(fine_grid_ptr, coarse_grid_ptr);
}

TEST(MeshMappingTest, Test2D)
{
  if (opensn::mpi_comm.size() != 1)
    return;

  const auto fine_grid_ptr =
    BuildOrthogonalMesh({{-1.0, -0.8, -0.75, -0.5, -0.25, 0.0, 0.5, 1.0},
                         {-1.0, -0.75, -0.5, -0.25, 0.0, 0.25, 0.5, 0.75, 1.0}});
  const auto coarse_grid_ptr =
    BuildOrthogonalMesh({{-1.0, -0.5, 0.0, 0.5, 1.0}, {-1.0, -0.5, 0.0, 0.5, 1.0}});
  TestMapping(fine_grid_ptr, coarse_grid_ptr);
}

TEST(MeshMappingTest, Test3D)
{
  if (opensn::mpi_comm.size() != 1)
    return;

  const auto fine_grid_ptr =
    BuildOrthogonalMesh({{-1.0, -0.75, -0.5, -0.25, 0.0, 0.5, 1.0},
                         {-1.0, -0.9, -0.75, -0.5, -0.25, 0.0, 0.25, 0.5, 0.75, 1.0},
                         {-1.0, -0.75, -0.5, -0.25, 0.0, 0.5, 1.0}});
  const auto coarse_grid_ptr =
    BuildOrthogonalMesh({{-1.0, 0.0, 1.0}, {-1.0, 0.0, 1.0}, {-1.0, -0.5, 0.0, 0.5, 1.0}});
  TestMapping(fine_grid_ptr, coarse_grid_ptr);
}

TEST(MeshMappingTest, ConcaveFaces)
{
  if (opensn::mpi_comm.size() != 1)
    return;

  // One cell mapped onto itself: the L-shaped polygon (0,0), (3,0), (3,1), (1,1), (1,3), (0,3)
  // extruded over z in [0, 1]. The centroid of its top and bottom faces, the average of their
  // vertices, is at (4/3, 4/3), outside the faces, so it cannot be the point that finds the
  // coarse face of a fine face.
  const auto file_name =
    std::filesystem::path(OPENSN_TEST_ROOT) / "assets" / "mesh" / "l_shaped_polygon.obj";
  const auto fine_grid_ptr = BuildExtrudedMeshFromFile(file_name, 1.0, 1);
  const auto coarse_grid_ptr = BuildExtrudedMeshFromFile(file_name, 1.0, 1);
  ASSERT_EQ(fine_grid_ptr->GetLocalCellCount(), std::size_t(1));
  ASSERT_EQ(coarse_grid_ptr->GetLocalCellCount(), std::size_t(1));
  const auto& fine_cell = fine_grid_ptr->GetLocalCell(0);
  const auto& coarse_cell = coarse_grid_ptr->GetLocalCell(0);
  ASSERT_EQ(fine_cell.faces.size(), std::size_t(8));
  ASSERT_EQ(coarse_cell.faces.size(), std::size_t(8));

  MeshMapping mesh_mapping;
  ASSERT_NO_THROW(mesh_mapping.Build(fine_grid_ptr, coarse_grid_ptr));

  const auto& fine_mapping = mesh_mapping.GetFineMapping(fine_cell);
  const auto& coarse_mapping = mesh_mapping.GetCoarseMapping(coarse_cell);
  EXPECT_EQ(fine_mapping.coarse_cell, &coarse_cell);
  ASSERT_EQ(coarse_mapping.fine_cells.size(), std::size_t(1));
  EXPECT_EQ(coarse_mapping.fine_cells[0], &fine_cell);

  // The two grids are the same, so every face maps to the face of the same index
  for (std::size_t face_i = 0; face_i < fine_cell.faces.size(); ++face_i)
  {
    EXPECT_EQ(fine_mapping.coarse_faces[face_i], face_i);
    ASSERT_EQ(coarse_mapping.fine_faces[face_i].size(), std::size_t(1));
    EXPECT_EQ(coarse_mapping.fine_faces[face_i][0], std::make_pair(&fine_cell, face_i));
  }
}

TEST(MeshMappingTest, GetCoarseMappingMissing)
{
  if (opensn::mpi_comm.size() != 1)
    return;

  const auto grid_ptr = BuildOrthogonalMesh({{-1.0, 1.0}});

  MeshMapping mesh_mapping;
  EXPECT_THROW(
    {
      try
      {
        for (const auto& cell : grid_ptr->GetLocalCells())
          mesh_mapping.GetCoarseMapping(*cell);
      }
      catch (const std::logic_error& e)
      {
        EXPECT_TRUE(std::string(e.what()).find("Coarse cell not found in mapping") !=
                    std::string::npos);
        throw;
      }
    },
    std::logic_error);
}

TEST(MeshMappingTest, GetFineMappingMissing)
{
  if (opensn::mpi_comm.size() != 1)
    return;

  const auto grid_ptr = BuildOrthogonalMesh({{-1.0, 1.0}});

  MeshMapping mesh_mapping;
  EXPECT_THROW(
    {
      try
      {
        for (const auto& cell : grid_ptr->GetLocalCells())
          mesh_mapping.GetFineMapping(*cell);
      }
      catch (const std::logic_error& e)
      {
        EXPECT_TRUE(std::string(e.what()).find("Fine cell not found in mapping") !=
                    std::string::npos);
        throw;
      }
    },
    std::logic_error);
}
