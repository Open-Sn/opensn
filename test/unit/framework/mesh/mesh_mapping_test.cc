#include "test/unit/opensn_unit_test.h"
#include "framework/mesh/mesh_mapping/mesh_mapping.h"

using namespace opensn;

class MeshMappingTest : public OpenSnUnitTest
{
};

void
TestMapping(const MeshContinuum& fine_grid, const MeshContinuum& coarse_grid)
{
  MeshMapping mesh_mapping;
  mesh_mapping.Build(fine_grid, coarse_grid);

  // Volumetric mapping; check all coarse cells against all fine cells
  for (const auto& coarse_cell : coarse_grid.local_cells)
  {
    const auto& coarse_mapping = mesh_mapping.GetCoarseMapping(coarse_cell);
    for (const auto& fine_cell : fine_grid.local_cells)
    {
      const auto in_coarse_mapping =
        std::find(coarse_mapping.fine_cells.begin(), coarse_mapping.fine_cells.end(), &fine_cell) !=
        coarse_mapping.fine_cells.end();
      const auto& fine_mapping = mesh_mapping.GetFineMapping(fine_cell);
      const auto in_fine_mapping = fine_mapping.coarse_cell == &coarse_cell;
      EXPECT_EQ(in_coarse_mapping, in_fine_mapping);

      const auto fine_within_coarse =
        coarse_grid.CheckPointInsideCell(coarse_cell, fine_cell.centroid);
      EXPECT_EQ(in_coarse_mapping, fine_within_coarse);
    }
  }

  // Surface mapping. Check all fine cell faces against all coarse cell faces
  for (const auto& fine_cell : fine_grid.local_cells)
  {
    const auto& fine_mapping = mesh_mapping.GetFineMapping(fine_cell);
    for (std::size_t fine_face_i = 0; fine_face_i < fine_cell.faces.size(); ++fine_face_i)
    {
      const auto& fine_face = fine_cell.faces[fine_face_i];
      bool found_mapping = false;
      for (const auto& coarse_cell : coarse_grid.local_cells)
      {
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
            coarse_grid.CheckPointInsideCellFace(coarse_cell, coarse_face_i, fine_face.centroid);
          const auto fine_within_coarse =
            coarse_grid.CheckPointInsideCell(coarse_cell, fine_cell.centroid);
          EXPECT_EQ(in_fine_mapping, fine_face_within_coarse && fine_within_coarse);
        }
      }
      if (!found_mapping)
        EXPECT_EQ(MeshMapping::invalid_face_index, fine_mapping.coarse_faces[fine_face_i]);
    }
  }
}

TEST_F(MeshMappingTest, Test1D)
{
  const auto fine_grid_ptr = BuildOrthogonalMesh({{-1.0, -0.8, -0.75, -0.5, -0.25, 0.0, 0.5, 1.0}});
  const auto coarse_grid_ptr = BuildOrthogonalMesh({{-1.0, -0.5, 0.0, 0.5, 1.0}});
  TestMapping(*fine_grid_ptr, *coarse_grid_ptr);
}

TEST_F(MeshMappingTest, Test2D)
{
  const auto fine_grid_ptr =
    BuildOrthogonalMesh({{-1.0, -0.8, -0.75, -0.5, -0.25, 0.0, 0.5, 1.0},
                         {-1.0, -0.75, -0.5, -0.25, 0.0, 0.25, 0.5, 0.75, 1.0}});
  const auto coarse_grid_ptr =
    BuildOrthogonalMesh({{-1.0, -0.5, 0.0, 0.5, 1.0}, {-1.0, -0.5, 0.0, 0.5, 1.0}});
  TestMapping(*fine_grid_ptr, *coarse_grid_ptr);
}

TEST_F(MeshMappingTest, Test3D)
{
  const auto fine_grid_ptr =
    BuildOrthogonalMesh({{-1.0, -0.75, -0.5, -0.25, 0.0, 0.5, 1.0},
                         {-1.0, -0.9, -0.75, -0.5, -0.25, 0.0, 0.25, 0.5, 0.75, 1.0},
                         {-1.0, -0.75, -0.5, -0.25, 0.0, 0.5, 1.0}});
  const auto coarse_grid_ptr =
    BuildOrthogonalMesh({{-1.0, 0.0, 1.0}, {-1.0, 0.0, 1.0}, {-1.0, -0.5, 0.0, 0.5, 1.0}});
  TestMapping(*fine_grid_ptr, *coarse_grid_ptr);
}

TEST_F(MeshMappingTest, GetCoarseMappingMissing)
{
  const auto grid_ptr = BuildOrthogonalMesh({{-1.0, 1.0}});

  MeshMapping mesh_mapping;
  EXPECT_THROW(
    {
      try
      {
        for (const auto& cell : grid_ptr->local_cells)
          mesh_mapping.GetCoarseMapping(cell);
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

TEST_F(MeshMappingTest, GetFineMappingMissing)
{
  const auto grid_ptr = BuildOrthogonalMesh({{-1.0, 1.0}});

  MeshMapping mesh_mapping;
  EXPECT_THROW(
    {
      try
      {
        for (const auto& cell : grid_ptr->local_cells)
          mesh_mapping.GetFineMapping(cell);
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
