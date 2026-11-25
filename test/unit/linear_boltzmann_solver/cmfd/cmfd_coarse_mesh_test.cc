#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/acceleration/cmfd_coarse_mesh.h"
#include "framework/mesh/logical_volume/logical_volume.h"
#include "framework/mesh/mesh/mesh.h"
#include "framework/mpi/mpi_utils.h"
#include "framework/runtime.h"
#include "test/unit/common/mesh_builders.h"
#include "gtest/gtest.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/acceleration/cmfd_vector_tools.h"
#include <algorithm>
#include <array>
#include <cmath>
#include <functional>
#include <limits>
#include <stdexcept>

using namespace opensn;

namespace
{

class PredicateLogicalVolume : public LogicalVolume
{
public:
  explicit PredicateLogicalVolume(std::function<bool(const Vector3&)> predicate)
    : predicate_(std::move(predicate))
  {
  }

  bool Inside(const Vector3& point) const override { return predicate_(point); }

private:
  std::function<bool(const Vector3&)> predicate_;
};

} // namespace

TEST(CMFDCoarseMesh, IdentityPreservesLocalCellGeometry)
{
  const unsigned int num_fine_cells =
    opensn::mpi_comm.size() > 2 ? static_cast<unsigned int>(opensn::mpi_comm.size()) : 2;
  auto grid = BuildLineMesh(static_cast<double>(num_fine_cells), num_fine_cells, 0.0);
  const auto coarse_mesh = CMFDCoarseMesh::BuildIdentity(*grid);

  ASSERT_EQ(coarse_mesh.NumLocalCells(), grid->GetLocalCellCount());
  ASSERT_EQ(coarse_mesh.LocalFineCellMemberships().size(), grid->GetLocalCellCount());

  for (const auto& fine_cell : grid->GetLocalCells())
  {
    ASSERT_TRUE(coarse_mesh.HasCoarseCell(fine_cell->global_id));
    EXPECT_EQ(coarse_mesh.MapFineCell(fine_cell->global_id), fine_cell->global_id);

    const auto& coarse_cell = coarse_mesh.LocalCell(fine_cell->local_id);
    EXPECT_EQ(coarse_cell.global_id, fine_cell->global_id);
    EXPECT_EQ(coarse_cell.local_id, fine_cell->local_id);
    EXPECT_EQ(coarse_cell.partition_id, fine_cell->partition_id);
    EXPECT_EQ(coarse_cell.block_id, fine_cell->block_id);
    EXPECT_DOUBLE_EQ(coarse_cell.volume, fine_cell->volume);
    EXPECT_EQ(coarse_cell.fine_cell_ids, std::vector<uint64_t>({fine_cell->global_id}));
    ASSERT_EQ(coarse_cell.faces.size(), fine_cell->faces.size());

    for (size_t f = 0; f < fine_cell->faces.size(); ++f)
    {
      const auto& coarse_face = coarse_cell.faces[f];
      const auto& fine_face = fine_cell->faces[f];
      EXPECT_EQ(coarse_face.has_neighbor, fine_face.has_neighbor);
      EXPECT_EQ(coarse_face.neighbor_id, fine_face.neighbor_id);
      EXPECT_DOUBLE_EQ(coarse_face.area, fine_face.area);
    }
  }

  for (const auto& membership : coarse_mesh.LocalFineCellMemberships())
  {
    const auto& fine_cell = grid->GetGlobalCell(membership.fine_cell_id);
    EXPECT_EQ(fine_cell.partition_id, opensn::mpi_comm.rank());
    EXPECT_EQ(membership.coarse_cell_id, fine_cell.global_id);
    EXPECT_EQ(membership.coarse_cell_partition_id, fine_cell.partition_id);
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
    (grid->GetLocalCellCount() + target_fine_cells_per_coarse_cell - 1) /
    target_fine_cells_per_coarse_cell;
  std::size_t expected_global_cells = 0;
  opensn::mpi_comm.all_reduce(
    expected_local_cells, expected_global_cells, mpi::op::sum<std::size_t>());

  ASSERT_EQ(coarse_mesh.NumGlobalCells(), expected_global_cells);
  ASSERT_EQ(coarse_mesh.NumLocalCells(), expected_local_cells);
  ASSERT_EQ(coarse_mesh.LocalFineCellMemberships().size(), grid->GetLocalCellCount());

  if (opensn::mpi_comm.size() > 1)
  {
    for (const auto& coarse_cell : coarse_mesh.LocalCells())
    {
      double expected_volume = 0.0;
      for (const auto fine_cell_id : coarse_cell.fine_cell_ids)
      {
        EXPECT_EQ(coarse_mesh.MapFineCell(fine_cell_id), coarse_cell.global_id);
        expected_volume += grid->GetGlobalCell(fine_cell_id).volume;
      }
      EXPECT_DOUBLE_EQ(coarse_cell.volume, expected_volume);
      EXPECT_FALSE(coarse_cell.faces.empty());
    }
    for (const auto& membership : coarse_mesh.LocalFineCellMemberships())
    {
      const auto& fine_cell = grid->GetGlobalCell(membership.fine_cell_id);
      EXPECT_EQ(fine_cell.partition_id, opensn::mpi_comm.rank());
      EXPECT_EQ(coarse_mesh.MapFineCell(fine_cell.global_id), membership.coarse_cell_id);
      EXPECT_EQ(membership.coarse_cell_partition_id, opensn::mpi_comm.rank());
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

  EXPECT_DOUBLE_EQ(first.volume, grid->GetGlobalCell(0).volume + grid->GetGlobalCell(1).volume);
  EXPECT_DOUBLE_EQ(second.volume, grid->GetGlobalCell(2).volume + grid->GetGlobalCell(3).volume);

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
    (grid->GetLocalCellCount() + target_fine_cells_per_coarse_cell - 1) /
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
        expected_area += grid->GetGlobalCell(fine_face.cell_id).faces[fine_face.face_index].area;

      EXPECT_DOUBLE_EQ(coarse_face.area, expected_area);
      found_merged_face = true;
    }
  }

  int global_found_merged_face = 0;
  opensn::mpi_comm.all_reduce(
    found_merged_face ? 1 : 0, global_found_merged_face, mpi::op::max<int>());
  EXPECT_TRUE(global_found_merged_face != 0);
}

TEST(CMFDCoarseMesh, MixedMaterialLocalAggregationCrossesBlockBoundaries)
{
  const auto num_fine_cells = 4U * static_cast<unsigned int>(opensn::mpi_comm.size());
  auto grid = BuildLineMesh(static_cast<double>(num_fine_cells), num_fine_cells, 0.0);
  for (auto& cell : grid->GetLocalCells())
    cell->block_id = static_cast<unsigned int>(cell->global_id % 2);
  for (const auto ghost_id : grid->GetGhostGlobalIDs())
  {
    auto& cell = grid->GetGlobalCell(ghost_id);
    cell.block_id = static_cast<unsigned int>(cell.global_id % 2);
  }

  const auto material_pure = CMFDCoarseMesh::BuildLocalAggregation(*grid, 4, false);
  const auto mixed = CMFDCoarseMesh::BuildLocalAggregation(*grid, 4, true);

  EXPECT_EQ(material_pure.NumLocalCells(), 4);
  ASSERT_EQ(mixed.NumLocalCells(), 1);
  EXPECT_EQ(mixed.LocalCell(0).fine_cell_ids.size(), 4);
  EXPECT_EQ(mixed.LocalFineCellMemberships().size(), 4);
}

TEST(CMFDCoarseMesh, GlobalAggregationBuildsConnectedCoarseCells)
{
  const std::size_t target_fine_cells_per_coarse_cell = 2;
  const unsigned int num_fine_cells =
    opensn::mpi_comm.size() > 2 ? 2 * static_cast<unsigned int>(opensn::mpi_comm.size()) : 4;
  auto grid = BuildLineMesh(static_cast<double>(num_fine_cells), num_fine_cells, 0.0);
  const auto coarse_mesh =
    CMFDCoarseMesh::BuildGlobalAggregation(*grid, target_fine_cells_per_coarse_cell);

  const std::size_t expected_global_cells =
    (grid->GetGlobalNumberOfCells() + target_fine_cells_per_coarse_cell - 1) /
    target_fine_cells_per_coarse_cell;
  ASSERT_EQ(coarse_mesh.NumGlobalCells(), expected_global_cells);
  ASSERT_EQ(coarse_mesh.LocalFineCellMemberships().size(), grid->GetLocalCellCount());

  for (const auto& membership : coarse_mesh.LocalFineCellMemberships())
  {
    const auto& fine_cell = grid->GetGlobalCell(membership.fine_cell_id);
    EXPECT_EQ(fine_cell.partition_id, opensn::mpi_comm.rank());
    EXPECT_TRUE(coarse_mesh.HasCoarseCell(fine_cell.global_id));
    EXPECT_EQ(coarse_mesh.MapFineCell(fine_cell.global_id), membership.coarse_cell_id);
  }

  std::size_t local_owned_fine_cells = 0;
  for (const auto& coarse_cell : coarse_mesh.LocalCells())
  {
    EXPECT_EQ(coarse_cell.partition_id, opensn::mpi_comm.rank());
    EXPECT_LE(coarse_cell.fine_cell_ids.size(), target_fine_cells_per_coarse_cell);
    EXPECT_FALSE(coarse_cell.fine_cell_ids.empty());
    EXPECT_GT(coarse_cell.volume, 0.0);
    local_owned_fine_cells += coarse_cell.fine_cell_ids.size();
  }

  std::size_t global_owned_fine_cells = 0;
  opensn::mpi_comm.all_reduce(
    local_owned_fine_cells, global_owned_fine_cells, mpi::op::sum<std::size_t>());
  EXPECT_EQ(global_owned_fine_cells, grid->GetGlobalNumberOfCells());
}

TEST(CMFDCoarseMesh, MixedMaterialGlobalAggregationCanMergeAcrossRanks)
{
  if (opensn::mpi_comm.size() < 2)
    GTEST_SKIP() << "Cross-rank aggregation requires at least two MPI ranks.";

  const auto num_fine_cells = 4U * static_cast<unsigned int>(opensn::mpi_comm.size());
  auto grid = BuildLineMesh(static_cast<double>(num_fine_cells), num_fine_cells, 0.0);
  for (auto& cell : grid->GetLocalCells())
    cell->block_id = static_cast<unsigned int>(cell->global_id % 2);
  for (const auto ghost_id : grid->GetGhostGlobalIDs())
  {
    auto& cell = grid->GetGlobalCell(ghost_id);
    cell.block_id = static_cast<unsigned int>(cell.global_id % 2);
  }

  const auto material_pure = CMFDCoarseMesh::BuildGlobalAggregation(*grid, 8, false);
  const auto mixed = CMFDCoarseMesh::BuildGlobalAggregation(*grid, 8, true);

  EXPECT_EQ(material_pure.NumGlobalCells(), num_fine_cells);
  EXPECT_EQ(mixed.NumGlobalCells(), (static_cast<std::size_t>(opensn::mpi_comm.size()) + 1) / 2);
}

TEST(CMFDCoarseMesh, LogicalVolumeOverlayMergesMaterialsAndIsDistributed)
{
  const auto num_fine_cells = 4U * static_cast<unsigned int>(opensn::mpi_comm.size());
  const auto length = static_cast<double>(num_fine_cells);
  auto grid = BuildLineMesh(length, num_fine_cells, 0.0);

  const auto SetBlock = [&](Cell& cell)
  { cell.block_id = cell.centroid.z < 0.25 * length ? 1U : 0U; };
  for (auto& cell : grid->GetLocalCells())
    SetBlock(*cell);
  for (const auto ghost_id : grid->GetGhostGlobalIDs())
    SetBlock(grid->GetGlobalCell(ghost_id));

  std::vector<std::shared_ptr<LogicalVolume>> regions;
  regions.push_back(std::make_shared<PredicateLogicalVolume>([length](const Vector3& point)
                                                             { return point.z < 0.5 * length; }));
  regions.push_back(std::make_shared<PredicateLogicalVolume>([length](const Vector3& point)
                                                             { return point.z >= 0.5 * length; }));

  const auto coarse_mesh = CMFDCoarseMesh::BuildLogicalVolumeOverlay(*grid, regions);
  ASSERT_EQ(coarse_mesh.NumGlobalCells(), 2);
  ASSERT_EQ(coarse_mesh.LocalFineCellMemberships().size(), grid->GetLocalCellCount());

  for (const auto& cell : grid->GetLocalCells())
  {
    const uint64_t expected_coarse_id = cell->centroid.z < 0.5 * length ? 0 : 1;
    EXPECT_EQ(coarse_mesh.MapFineCell(cell->global_id), expected_coarse_id);
  }

  const std::array<double, 2> expected_volumes = {0.5 * length, 0.5 * length};
  for (uint64_t coarse_id = 0; coarse_id < expected_volumes.size(); ++coarse_id)
  {
    const double local_volume = coarse_mesh.HasLocalCoarseCell(coarse_id)
                                  ? coarse_mesh.LocalCellFromGlobalID(coarse_id).volume
                                  : 0.0;
    double global_volume = 0.0;
    opensn::mpi_comm.all_reduce(local_volume, global_volume, mpi::op::sum<double>());
    EXPECT_DOUBLE_EQ(global_volume, expected_volumes[coarse_id]);
  }

  std::size_t local_internal_faces = 0;
  for (const auto& coarse_cell : coarse_mesh.LocalCells())
    local_internal_faces +=
      static_cast<std::size_t>(std::count_if(coarse_cell.faces.begin(),
                                             coarse_cell.faces.end(),
                                             [](const auto& face) { return face.has_neighbor; }));
  std::size_t global_internal_faces = 0;
  opensn::mpi_comm.all_reduce(
    local_internal_faces, global_internal_faces, mpi::op::sum<std::size_t>());
  EXPECT_EQ(global_internal_faces, 2);
}

TEST(CMFDCoarseMesh, LogicalVolumeOverlayRejectsHolesAndOverlaps)
{
  auto grid = BuildLineMesh(4.0, 4, 0.0);
  const auto all = std::make_shared<PredicateLogicalVolume>([](const Vector3&) { return true; });
  const auto left =
    std::make_shared<PredicateLogicalVolume>([](const Vector3& point) { return point.z < 2.0; });

  EXPECT_THROW(CMFDCoarseMesh::BuildLogicalVolumeOverlay(*grid, {all, all}), std::logic_error);
  EXPECT_THROW(CMFDCoarseMesh::BuildLogicalVolumeOverlay(*grid, {left}), std::logic_error);
}

TEST(CMFDCoarseMesh, LogicalVolumeOverlaySplitsDisconnectedRegion)
{
  auto grid = BuildLineMesh(4.0, 4, 0.0);
  const auto ends = std::make_shared<PredicateLogicalVolume>(
    [](const Vector3& point) { return point.z < 1.0 or point.z >= 3.0; });
  const auto middle = std::make_shared<PredicateLogicalVolume>(
    [](const Vector3& point) { return point.z >= 1.0 and point.z < 3.0; });

  const auto coarse_mesh = CMFDCoarseMesh::BuildLogicalVolumeOverlay(*grid, {ends, middle});
  EXPECT_EQ(coarse_mesh.NumGlobalCells(), 3);
  for (const auto& cell : grid->GetLocalCells())
  {
    const auto coarse_id = coarse_mesh.MapFineCell(cell->global_id);
    if (cell->centroid.z < 1.0)
      EXPECT_EQ(coarse_id, 0);
    else if (cell->centroid.z >= 3.0)
      EXPECT_EQ(coarse_id, 1);
    else
      EXPECT_EQ(coarse_id, 2);
  }
}

TEST(CMFDCoarseMesh, CartesianOverlayUsesHalfOpenIntervalsAndValidatesCoverage)
{
  auto grid = BuildLineMesh(4.0, 4, 0.0);
  const auto coarse_mesh = CMFDCoarseMesh::BuildCartesianOverlay(*grid, {}, {}, {0.0, 1.5, 4.0});
  EXPECT_EQ(coarse_mesh.NumGlobalCells(), 2);
  for (const auto& cell : grid->GetLocalCells())
  {
    const uint64_t expected_id = cell->centroid.z < 1.5 ? 0 : 1;
    EXPECT_EQ(coarse_mesh.MapFineCell(cell->global_id), expected_id);
  }

  EXPECT_THROW(CMFDCoarseMesh::BuildCartesianOverlay(*grid, {}, {}, {}), std::invalid_argument);
  EXPECT_THROW(CMFDCoarseMesh::BuildCartesianOverlay(*grid, {}, {}, {0.0, 1.0, 1.0}),
               std::invalid_argument);
  EXPECT_THROW(CMFDCoarseMesh::BuildCartesianOverlay(*grid, {}, {}, {0.0, 3.0}), std::logic_error);
}

TEST(CMFDCoarseMesh, CartesianOverlayReportsConformingInterfaces)
{
  // A 4x4 mesh of unit squares split into 2x2 overlay regions along x = 2 and y = 2. Each plane
  // carries four unit-length fine faces, and every one of them lies on the plane.
  auto grid = BuildSquareMesh(4.0, 4, 0.0);
  const auto coarse_mesh =
    CMFDCoarseMesh::BuildCartesianOverlay(*grid, {0.0, 2.0, 4.0}, {0.0, 2.0, 4.0}, {});
  EXPECT_EQ(coarse_mesh.NumGlobalCells(), 4);
  const auto& overlay_diagnostics = coarse_mesh.OverlayDiagnostics();
  ASSERT_TRUE(overlay_diagnostics.has_value());
  const auto diagnostics = overlay_diagnostics.value_or(CMFDOverlayDiagnostics{});
  EXPECT_EQ(diagnostics.num_interface_faces, 8);
  EXPECT_EQ(diagnostics.num_nonconforming_faces, 0);
  EXPECT_DOUBLE_EQ(diagnostics.interface_area, 8.0);
  EXPECT_DOUBLE_EQ(diagnostics.nonconforming_area, 0.0);
  EXPECT_EQ(diagnostics.num_thin_intervals, 0);
  EXPECT_DOUBLE_EQ(diagnostics.min_interval_width_ratio, 1.0);
}

TEST(CMFDCoarseMesh, CartesianOverlayReportsNonconformingInterfaces)
{
  // The plane x = 1.5 passes through the second column of cells. Cells are assigned by centroid,
  // so the region interface follows the fine faces at x = 1, which do not lie on the plane.
  auto grid = BuildSquareMesh(4.0, 4, 0.0);
  const auto coarse_mesh = CMFDCoarseMesh::BuildCartesianOverlay(*grid, {0.0, 1.5, 4.0}, {}, {});
  EXPECT_EQ(coarse_mesh.NumGlobalCells(), 2);
  const auto& overlay_diagnostics = coarse_mesh.OverlayDiagnostics();
  ASSERT_TRUE(overlay_diagnostics.has_value());
  const auto diagnostics = overlay_diagnostics.value_or(CMFDOverlayDiagnostics{});
  EXPECT_EQ(diagnostics.num_interface_faces, 4);
  EXPECT_EQ(diagnostics.num_nonconforming_faces, 4);
  EXPECT_DOUBLE_EQ(diagnostics.nonconforming_area, 4.0);
  EXPECT_EQ(diagnostics.num_thin_intervals, 0);
}

TEST(CMFDCoarseMesh, CartesianOverlayReportsThinOccupiedIntervals)
{
  // The interval [0.45, 0.55) holds the centroid of the first cell and is 0.1 wide, against an
  // upper-median occupied width of 3.45. The empty interval [0, 0.45) is not counted.
  auto grid = BuildLineMesh(4.0, 4, 0.0);
  const auto coarse_mesh =
    CMFDCoarseMesh::BuildCartesianOverlay(*grid, {}, {}, {0.0, 0.45, 0.55, 4.0});
  EXPECT_EQ(coarse_mesh.NumGlobalCells(), 2);
  const auto& overlay_diagnostics = coarse_mesh.OverlayDiagnostics();
  ASSERT_TRUE(overlay_diagnostics.has_value());
  const auto diagnostics = overlay_diagnostics.value_or(CMFDOverlayDiagnostics{});
  EXPECT_EQ(diagnostics.num_thin_intervals, 1);
  EXPECT_NEAR(diagnostics.min_interval_width_ratio, 0.1 / 3.45, 1.0e-12);
  EXPECT_EQ(diagnostics.num_interface_faces, 1);
  EXPECT_EQ(diagnostics.num_nonconforming_faces, 1);
}

TEST(CMFDCoarseMesh, LogicalVolumeOverlayHasNoCartesianDiagnostics)
{
  auto grid = BuildLineMesh(4.0, 4, 0.0);
  const auto all = std::make_shared<PredicateLogicalVolume>([](const Vector3&) { return true; });
  const auto coarse_mesh = CMFDCoarseMesh::BuildLogicalVolumeOverlay(*grid, {all});
  EXPECT_EQ(coarse_mesh.NumGlobalCells(), 1);
  EXPECT_FALSE(coarse_mesh.OverlayDiagnostics().has_value());
}

TEST(CMFDBalance, RelativeResidualIsIndependentOfFluxScale)
{
  for (const double scale : {1.0, 1.0e-200, 1.0e200})
  {
    std::vector<std::pair<double, double>> values;
    // Also exercise ranks with no coarse cells.
    if (mpi_comm.rank() == 0)
      values = {{3.0 * scale, 6.0 * scale}, {4.0 * scale, 8.0 * scale}};
    EXPECT_DOUBLE_EQ(CMFDRelativeBalanceResidual(values), 0.5);
  }
}

TEST(CMFDBalance, DegenerateAndNonfiniteBalancesCannotConverge)
{
  EXPECT_TRUE(std::isinf(CMFDRelativeBalanceResidual({{0.0, 0.0}})));
  EXPECT_TRUE(std::isinf(CMFDRelativeBalanceResidual({{1.0, 0.0}})));
  std::vector<std::pair<double, double>> values = {{0.0, 1.0}};
  if (mpi_comm.rank() == 0)
    values[0].first = std::numeric_limits<double>::quiet_NaN();
  EXPECT_TRUE(std::isinf(CMFDRelativeBalanceResidual(values)));
  EXPECT_DOUBLE_EQ(CMFDRelativeBalanceResidual({{0.0, 1.0e-200}}), 0.0);
}

TEST(CMFDCorrection, SignChangingRatiosAreCountedAndExistingSignsAreKept)
{
  const double nan = std::numeric_limits<double>::quiet_NaN();
  const double inf = std::numeric_limits<double>::infinity();
  // Positive uncorrected values that the correction only rescales.
  EXPECT_EQ(CMFDCountSignChangingRatios({1.0, 3.0}, {5.0, 1.0e-300}), 0U);
  // A negative uncorrected value is not rescaled, even by a positive ratio, which could amplify
  // it.
  EXPECT_EQ(CMFDCountSignChangingRatios({-2.0, -1.0e-3}, {-0.5, -10.0}), 2U);
  // Flipped signs and non-finite corrected values change or destroy the fine flux.
  EXPECT_EQ(CMFDCountSignChangingRatios({1.0, -2.0, 4.0, 5.0}, {-1.0, 2.0, nan, inf}), 4U);
  // A zero corrected value sets the fine flux to zero without a sign change, as the coarse
  // solution of a decoupled source-free region does. Negative zero is zero.
  EXPECT_EQ(CMFDCountSignChangingRatios({3.0, 3.0}, {0.0, -0.0}), 0U);
  // A negligible uncorrected value is not rescaled, so its corrected value is irrelevant.
  EXPECT_EQ(CMFDCountSignChangingRatios({0.0, 1.0e-15, -1.0e-15}, {-1.0, -1.0, 1.0}), 0U);
  // A non-finite uncorrected value has no admissible ratio.
  EXPECT_EQ(CMFDCountSignChangingRatios({nan}, {1.0}), 1U);
  EXPECT_EQ(CMFDCountSignChangingRatios({}, {}), 0U);
  EXPECT_THROW(CMFDCountSignChangingRatios({1.0}, {1.0, 2.0}), std::invalid_argument);
}
