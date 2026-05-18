#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/acceleration/cmfd_coarse_mesh.h"
#include "framework/mpi/mpi_utils.h"
#include "framework/runtime.h"
#include "test/unit/common/mesh_builders.h"
#include "gtest/gtest.h"

using namespace opensn;

TEST(CMFDCoarseMesh, IdentityPreservesLocalCellGeometry)
{
  const unsigned int num_fine_cells =
    opensn::mpi_comm.size() > 2 ? static_cast<unsigned int>(opensn::mpi_comm.size()) : 2;
  auto grid = BuildLineMesh(static_cast<double>(num_fine_cells), num_fine_cells, 0.0);
  const auto coarse_mesh = CMFDCoarseMesh::BuildIdentity(*grid);

  ASSERT_EQ(coarse_mesh.NumLocalCells(), grid->local_cells.size());

  for (const auto& fine_cell : grid->local_cells)
  {
    ASSERT_TRUE(coarse_mesh.HasCoarseCell(fine_cell.global_id));
    EXPECT_EQ(coarse_mesh.MapFineCell(fine_cell.global_id), fine_cell.global_id);

    const auto& coarse_cell = coarse_mesh.LocalCell(fine_cell.local_id);
    EXPECT_EQ(coarse_cell.global_id, fine_cell.global_id);
    EXPECT_EQ(coarse_cell.local_id, fine_cell.local_id);
    EXPECT_EQ(coarse_cell.partition_id, fine_cell.partition_id);
    EXPECT_EQ(coarse_cell.block_id, fine_cell.block_id);
    EXPECT_DOUBLE_EQ(coarse_cell.volume, fine_cell.volume);
    EXPECT_EQ(coarse_cell.fine_cell_ids, std::vector<uint64_t>({fine_cell.global_id}));
    ASSERT_EQ(coarse_cell.faces.size(), fine_cell.faces.size());

    for (size_t f = 0; f < fine_cell.faces.size(); ++f)
    {
      const auto& coarse_face = coarse_cell.faces[f];
      const auto& fine_face = fine_cell.faces[f];
      EXPECT_EQ(coarse_face.has_neighbor, fine_face.has_neighbor);
      EXPECT_EQ(coarse_face.neighbor_id, fine_face.neighbor_id);
      EXPECT_DOUBLE_EQ(coarse_face.area, fine_face.area);
    }
  }
}

TEST(CMFDCoarseMesh, LocalAggregationBuildsConnectedCoarseCells)
{
  const std::size_t target_fine_cells_per_coarse_cell = 2;
  const unsigned int num_fine_cells =
    opensn::mpi_comm.size() > 2 ? 2 * static_cast<unsigned int>(opensn::mpi_comm.size()) : 4;
  auto grid = BuildLineMesh(static_cast<double>(num_fine_cells), num_fine_cells, 0.0);
  const auto coarse_mesh =
    CMFDCoarseMesh::BuildLocalAggregation(*grid, target_fine_cells_per_coarse_cell);

  const std::size_t expected_local_cells =
    (grid->local_cells.size() + target_fine_cells_per_coarse_cell - 1) /
    target_fine_cells_per_coarse_cell;
  std::size_t expected_global_cells = 0;
  opensn::mpi_comm.all_reduce(
    expected_local_cells, expected_global_cells, mpi::op::sum<std::size_t>());

  ASSERT_EQ(coarse_mesh.NumGlobalCells(), expected_global_cells);
  ASSERT_EQ(coarse_mesh.NumLocalCells(), expected_local_cells);

  if (opensn::mpi_comm.size() > 1)
  {
    for (const auto& coarse_cell : coarse_mesh.LocalCells())
    {
      double expected_volume = 0.0;
      for (const auto fine_cell_id : coarse_cell.fine_cell_ids)
      {
        EXPECT_EQ(coarse_mesh.MapFineCell(fine_cell_id), coarse_cell.global_id);
        expected_volume += grid->cells[fine_cell_id].volume;
      }
      EXPECT_DOUBLE_EQ(coarse_cell.volume, expected_volume);
      EXPECT_FALSE(coarse_cell.faces.empty());
    }
    return;
  }

  ASSERT_EQ(coarse_mesh.NumLocalCells(), 2);
  const auto& first = coarse_mesh.LocalCell(0);
  const auto& second = coarse_mesh.LocalCell(1);

  EXPECT_EQ(first.fine_cell_ids, std::vector<uint64_t>({0, 1}));
  EXPECT_EQ(second.fine_cell_ids, std::vector<uint64_t>({2, 3}));
  EXPECT_EQ(coarse_mesh.MapFineCell(0), first.global_id);
  EXPECT_EQ(coarse_mesh.MapFineCell(1), first.global_id);
  EXPECT_EQ(coarse_mesh.MapFineCell(2), second.global_id);
  EXPECT_EQ(coarse_mesh.MapFineCell(3), second.global_id);

  EXPECT_DOUBLE_EQ(first.volume, grid->cells[0].volume + grid->cells[1].volume);
  EXPECT_DOUBLE_EQ(second.volume, grid->cells[2].volume + grid->cells[3].volume);

  ASSERT_EQ(first.faces.size(), 2);
  ASSERT_EQ(second.faces.size(), 2);

  bool first_has_boundary = false;
  bool first_has_second = false;
  for (const auto& face : first.faces)
  {
    first_has_boundary = first_has_boundary or not face.has_neighbor;
    first_has_second =
      first_has_second or (face.has_neighbor and face.neighbor_id == second.global_id);
  }

  bool second_has_boundary = false;
  bool second_has_first = false;
  for (const auto& face : second.faces)
  {
    second_has_boundary = second_has_boundary or not face.has_neighbor;
    second_has_first =
      second_has_first or (face.has_neighbor and face.neighbor_id == first.global_id);
  }

  EXPECT_TRUE(first_has_boundary);
  EXPECT_TRUE(first_has_second);
  EXPECT_TRUE(second_has_boundary);
  EXPECT_TRUE(second_has_first);
}

TEST(CMFDCoarseMesh, LocalAggregationMergesFineFacesOnSameCoarseInterface)
{
  const std::size_t target_fine_cells_per_coarse_cell = 2;
  const unsigned int num_cells_per_dimension = opensn::mpi_comm.size() > 1 ? 4 : 2;
  auto grid =
    BuildSquareMesh(static_cast<double>(num_cells_per_dimension), num_cells_per_dimension, 0.0);
  const auto coarse_mesh =
    CMFDCoarseMesh::BuildLocalAggregation(*grid, target_fine_cells_per_coarse_cell);

  const std::size_t expected_local_cells =
    (grid->local_cells.size() + target_fine_cells_per_coarse_cell - 1) /
    target_fine_cells_per_coarse_cell;
  std::size_t expected_global_cells = 0;
  opensn::mpi_comm.all_reduce(
    expected_local_cells, expected_global_cells, mpi::op::sum<std::size_t>());

  ASSERT_EQ(coarse_mesh.NumGlobalCells(), expected_global_cells);
  ASSERT_EQ(coarse_mesh.NumLocalCells(), expected_local_cells);

  bool found_merged_face = false;
  for (const auto& coarse_cell : coarse_mesh.LocalCells())
  {
    for (const auto& coarse_face : coarse_cell.faces)
    {
      if (coarse_face.fine_faces.size() < 2)
        continue;

      double expected_area = 0.0;
      for (const auto& fine_face : coarse_face.fine_faces)
        expected_area += grid->cells[fine_face.cell_id].faces[fine_face.face_index].area;

      EXPECT_DOUBLE_EQ(coarse_face.area, expected_area);
      found_merged_face = true;
    }
  }

  int global_found_merged_face = 0;
  opensn::mpi_comm.all_reduce(
    found_merged_face ? 1 : 0, global_found_merged_face, mpi::op::max<int>());
  EXPECT_TRUE(global_found_merged_face != 0);
}
