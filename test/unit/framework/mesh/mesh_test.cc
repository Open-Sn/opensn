#include "gtest/gtest.h"
#include "test/unit/common/mesh_builders.h"
#include "framework/mesh/mesh/mesh.h"
#include "framework/mesh/logical_volume/rpp_logical_volume.h"
#include "framework/parameters/parameter_block.h"
#include "framework/runtime.h"
#include <cmath>

using namespace opensn;

/// Helper for the PointInsideCellXD tests
void
TestPointInsideCell(const std::shared_ptr<Mesh> grid)
{
  // Centroid is contained within cell whose centroid it is
  for (const auto& cell : grid->GetLocalCells())
    for (const auto& other_cell : grid->GetLocalCells())
    {
      const auto same_cell = cell->global_id == other_cell->global_id;
      const auto within = grid->CheckPointInsideCell(*other_cell, cell->centroid);
      EXPECT_EQ(same_cell, within);
    }

  // Vertices are contained within cells
  for (const auto& cell : grid->GetLocalCells())
    for (const auto vi : cell->vertex_ids)
      for (const auto& other_cell : grid->GetLocalCells())
      {
        const auto has_vertex =
          std::find(other_cell->vertex_ids.begin(), other_cell->vertex_ids.end(), vi) !=
          other_cell->vertex_ids.end();
        const auto within = grid->CheckPointInsideCell(*other_cell, grid->GlobalVertex(vi));
        EXPECT_EQ(has_vertex, within);
      }

  // Face centroids are contained within cells (including neighbors)
  for (const auto& cell : grid->GetLocalCells())
    for (const auto& face : cell->faces)
      for (const auto& other_cell : grid->GetLocalCells())
      {
        const auto same_cell_or_neighbor =
          (cell->global_id == other_cell->global_id ||
           (face.has_neighbor && face.neighbor_id == other_cell->global_id));
        const auto within = grid->CheckPointInsideCell(*other_cell, face.centroid);
        EXPECT_EQ(same_cell_or_neighbor, within);
      }

  if (grid->GetDimension() > 1)
  {
    // Face edge centroids contained
    for (const auto& cell : grid->GetLocalCells())
      for (const auto& face : cell->faces)
        for (size_t side = 0; side < face.vertex_ids.size(); ++side)
        {
          const size_t sp1 = (side < (face.vertex_ids.size() - 1)) ? side + 1 : 0;
          const auto& v0 = grid->GlobalVertex(face.vertex_ids[side]);
          const auto& v1 = grid->GlobalVertex(face.vertex_ids[sp1]);
          const auto c = (v0 + v1) / 2.0;
          EXPECT_TRUE(grid->CheckPointInsideCell(*cell, c));
        }
  }

  if (grid->GetDimension() == 3)
  {
    // Tetrahedral centroids contained
    for (const auto& cell : grid->GetLocalCells())
    {
      for (const auto& face : cell->faces)
        for (size_t side = 0; side < face.vertex_ids.size(); ++side)
        {
          const size_t sp1 = (side < (face.vertex_ids.size() - 1)) ? side + 1 : 0;
          const auto& v0 = grid->GlobalVertex(face.vertex_ids[side]);
          const auto& v1 = face.centroid;
          const auto& v2 = grid->GlobalVertex(face.vertex_ids[sp1]);
          const auto& v3 = cell->centroid;
          const auto c = (v0 + v1 + v2 + v3) / 4.0;
          for (const auto& other_cell : grid->GetLocalCells())
          {
            const auto same_cell = cell->global_id == other_cell->global_id;
            const auto within = grid->CheckPointInsideCell(*other_cell, c);
            EXPECT_EQ(same_cell, within);
          }
        }
    }

    // Tetrahedral face centroids contained
    for (const auto& cell : grid->GetLocalCells())
      for (const auto& face : cell->faces)
        for (size_t side = 0; side < face.vertex_ids.size(); ++side)
        {
          const auto tet_face_vertices = grid->GetTetrahedralFaceVertices(*cell, face, side);
          for (const auto& v : tet_face_vertices)
            for (const auto& other_cell : grid->GetLocalCells())
            {
              const auto same_cell_or_neighbor =
                (cell->global_id == other_cell->global_id ||
                 (face.has_neighbor && face.neighbor_id == other_cell->global_id));
              const auto within = grid->CheckPointInsideCell(*other_cell, face.centroid);
              EXPECT_EQ(same_cell_or_neighbor, within);
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
      for (std::size_t face_i = 0; face_i < cell->faces.size(); ++face_i)
      {
        const auto& face = cell->faces[face_i];
        const auto has_vertex =
          std::find(face.vertex_ids.begin(), face.vertex_ids.end(), vi) != face.vertex_ids.end();
        const auto within = grid->CheckPointInsideCellFace(*cell, face_i, grid->GlobalVertex(vi));
        EXPECT_EQ(has_vertex, within);
      }

  // Cell centroids not contained within any faces
  for (const auto& cell : grid->GetLocalCells())
    for (const auto& other_cell : grid->GetLocalCells())
    {
      for (std::size_t other_face_i = 0; other_face_i < other_cell->faces.size(); ++other_face_i)
        EXPECT_FALSE(grid->CheckPointInsideCellFace(*other_cell, other_face_i, cell->centroid));
    }

  if (grid->GetDimension() > 1)
  {
    // Face centroids contained within faces
    for (const auto& cell : grid->GetLocalCells())
      for (const auto& face : cell->faces)
        for (const auto& other_cell : grid->GetLocalCells())
          for (std::size_t other_face_i = 0; other_face_i < other_cell->faces.size();
               ++other_face_i)
          {
            const auto& other_face = other_cell->faces[other_face_i];
            const auto same_face = other_face.centroid.AbsoluteEquals(face.centroid);
            const auto within =
              grid->CheckPointInsideCellFace(*other_cell, other_face_i, face.centroid);
            EXPECT_EQ(same_face, within);
          }
  }

  if (grid->GetDimension() > 2)
  {
    // Face edge centroids contained within faces that have the vertices that contain the edge
    for (const auto& cell : grid->GetLocalCells())
      for (const auto& face : cell->faces)
        for (std::size_t i = 0; i < face.vertex_ids.size(); ++i)
        {
          const auto vi1 = face.vertex_ids[i];
          const auto vi2 =
            (i == face.vertex_ids.size() - 1) ? face.vertex_ids[0] : face.vertex_ids[i + 1];
          const auto edge_centroid = (grid->GlobalVertex(vi1) + grid->GlobalVertex(vi2)) / 2;
          for (const auto& other_cell : grid->GetLocalCells())
            for (std::size_t other_face_i = 0; other_face_i < other_cell->faces.size();
                 ++other_face_i)
            {
              const auto& other_face = other_cell->faces[other_face_i];
              const auto has_edges =
                (std::find(other_face.vertex_ids.begin(), other_face.vertex_ids.end(), vi1) !=
                 other_face.vertex_ids.end()) &&
                (std::find(other_face.vertex_ids.begin(), other_face.vertex_ids.end(), vi2) !=
                 other_face.vertex_ids.end());
              const auto within =
                grid->CheckPointInsideCellFace(*other_cell, other_face_i, edge_centroid);
              EXPECT_EQ(has_edges, within);
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

namespace
{

/// RPP logical volume that is infinite in x and y and spans [zmin, zmax] in z.
std::shared_ptr<RPPLogicalVolume>
MakeZSlabVolume(double zmin, double zmax)
{
  auto params = RPPLogicalVolume::GetInputParameters();
  ParameterBlock block;
  block.AddParameter("infx", true);
  block.AddParameter("infy", true);
  block.AddParameter("zmin", zmin);
  block.AddParameter("zmax", zmax);
  params.AssignParameters(block);
  return std::make_shared<RPPLogicalVolume>(params);
}

/// Global number of boundary faces with their centroid at z that are named `name`.
int
CountNamedBoundaryFaces(const Mesh& grid, const std::string& name, double z)
{
  int local_count = 0;
  const auto& id_map = grid.GetBoundaryIDMap();
  for (const auto& cell : grid.GetLocalCells())
    for (const auto& face : cell->faces)
    {
      if (face.has_neighbor or std::fabs(face.centroid.z - z) > 1.0e-12)
        continue;
      const auto it = id_map.find(face.neighbor_id);
      if (it != id_map.end() and it->second == name)
        ++local_count;
    }
  int global_count = 0;
  mpi_comm.all_reduce(local_count, global_count, mpi::op::sum<int>());
  return global_count;
}

/// Checks that each ghost cell has `outside_id` if its centroid has z > 1 and `inside_id`
/// otherwise. Returns the global number of ghost cells with z > 1.
int
CheckGhostBlockIDs(const Mesh& grid, unsigned int outside_id, unsigned int inside_id)
{
  int local_count = 0;
  for (const auto ghost_id : grid.GetGhostGlobalIDs())
  {
    const auto& cell = grid.GetGlobalCell(ghost_id);
    const bool outside = cell.centroid.z > 1.0;
    EXPECT_EQ(cell.block_id, outside ? outside_id : inside_id);
    local_count += outside ? 1 : 0;
  }
  int global_count = 0;
  mpi_comm.all_reduce(local_count, global_count, mpi::op::sum<int>());
  return global_count;
}

} // namespace

TEST(MeshTest, SetBlockIDFromLogicalVolumeSense)
{
  // Eight cells of length 0.5 on [0, 4]; the volume holds the centroids of the first two.
  const auto grid = BuildLineMesh(4.0, 8, 0.0);
  grid->SetUniformBlockID(0);
  const auto volume = MakeZSlabVolume(0.0, 1.0);

  // sense = false selects the cells whose centroids are outside the volume, ghost cells too.
  grid->SetBlockIDFromLogicalVolume(*volume, 7, false);
  auto block_volumes = grid->ComputeVolumePerBlockID();
  EXPECT_EQ(block_volumes.size(), 2U);
  EXPECT_NEAR(block_volumes[0], 1.0, 1.0e-12);
  EXPECT_NEAR(block_volumes[7], 3.0, 1.0e-12);
  // Ghost cells outside the volume exist exactly when there are two or more ranks.
  EXPECT_EQ(CheckGhostBlockIDs(*grid, 7, 0) > 0, mpi_comm.size() > 1);

  // sense = true selects the cells whose centroids are inside the volume.
  grid->SetBlockIDFromLogicalVolume(*volume, 3, true);
  block_volumes = grid->ComputeVolumePerBlockID();
  EXPECT_EQ(block_volumes.size(), 2U);
  EXPECT_NEAR(block_volumes[3], 1.0, 1.0e-12);
  EXPECT_NEAR(block_volumes[7], 3.0, 1.0e-12);
  EXPECT_EQ(CheckGhostBlockIDs(*grid, 7, 3) > 0, mpi_comm.size() > 1);
}

TEST(MeshTest, SetBoundaryIDFromLogicalVolumeSense)
{
  // The boundary faces of the line mesh are at z = 0 and z = 4; the volume holds z = 0 only.
  const auto grid = BuildLineMesh(4.0, 8, 0.0);
  const auto volume = MakeZSlabVolume(-0.5, 0.5);

  // sense = false selects the boundary faces whose centroids are outside the volume.
  grid->SetBoundaryIDFromLogicalVolume(*volume, "outside_face", false);
  EXPECT_EQ(CountNamedBoundaryFaces(*grid, "outside_face", 0.0), 0);
  EXPECT_EQ(CountNamedBoundaryFaces(*grid, "outside_face", 4.0), 1);

  // sense = true selects the boundary faces whose centroids are inside the volume.
  grid->SetBoundaryIDFromLogicalVolume(*volume, "inside_face", true);
  EXPECT_EQ(CountNamedBoundaryFaces(*grid, "inside_face", 0.0), 1);
  EXPECT_EQ(CountNamedBoundaryFaces(*grid, "outside_face", 4.0), 1);
}

TEST(MeshTest, SetBoundaryIDFromInteriorLogicalVolumeSense)
{
  // The volume holds no boundary face, so sense = false selects both and sense = true none.
  const auto grid = BuildLineMesh(4.0, 8, 0.0);
  const auto volume = MakeZSlabVolume(1.0, 3.0);

  grid->SetBoundaryIDFromLogicalVolume(*volume, "outer", false);
  EXPECT_EQ(grid->GetBoundaryNameMap().count("outer"), 1U);
  EXPECT_EQ(CountNamedBoundaryFaces(*grid, "outer", 0.0), 1);
  EXPECT_EQ(CountNamedBoundaryFaces(*grid, "outer", 4.0), 1);

  grid->SetBoundaryIDFromLogicalVolume(*volume, "inner", true);
  EXPECT_EQ(grid->GetBoundaryNameMap().count("inner"), 0U);
  EXPECT_EQ(CountNamedBoundaryFaces(*grid, "outer", 0.0), 1);
  EXPECT_EQ(CountNamedBoundaryFaces(*grid, "outer", 4.0), 1);
}
