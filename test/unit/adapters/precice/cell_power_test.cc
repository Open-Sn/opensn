// SPDX-FileCopyrightText: 2026 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#ifdef OPENSN_WITH_PRECICE

#include "adapters/precice/cell_power.h"
#include "framework/field_functions/field_function_grid_based.h"
#include "framework/math/spatial_discretization/finite_element/piecewise_linear/piecewise_linear_discontinuous.h"
#include "framework/runtime.h"
#include "test/unit/common/mesh_builders.h"
#include "gtest/gtest.h"
#include <cmath>
#include <functional>
#include <limits>

using namespace opensn;

namespace
{
// Interpolate `function` at the PWLD nodes. Each test field lies in the PWLD space, so its cell
// integrals are exact up to floating-point rounding.
std::shared_ptr<FieldFunctionGridBased>
MakeNodalField(const std::shared_ptr<Mesh>& grid,
               const std::function<double(const Vector3&)>& function)
{
  std::shared_ptr<SpatialDiscretization> sdm = PieceWiseLinearDiscontinuous::New(grid);
  auto field =
    std::make_shared<FieldFunctionGridBased>("deposition", sdm, Unknown(UnknownType::SCALAR));
  const auto& unknowns = field->GetUnknownManager();
  std::vector<double> values(sdm->GetNumLocalDOFs(unknowns));
  for (const auto& cell : grid->GetLocalCells())
  {
    const auto& nodes = sdm->GetCellNodeLocations(*cell);
    for (std::size_t node = 0; node < nodes.size(); ++node)
      values[sdm->MapDOFLocal(*cell, node, unknowns, 0, 0)] = function(nodes[node]);
  }
  field->UpdateFieldVector(values);
  return field;
}

double
GlobalPower(const coupling::CellPowerData& cell_data)
{
  double global_power = 0.0;
  mpi_comm.all_reduce(cell_data.GetLocalPower(), global_power, mpi::op::sum<double>());
  return global_power;
}
} // namespace

TEST(PreciceCellPower, IntegratesOwnedLinearField1D)
{
  auto grid = BuildLineMesh(4.0, 4, 0.0);
  for (const auto& cell : grid->GetLocalCells())
    cell->block_id = 7;
  const auto field = MakeNodalField(grid, [](const Vector3& x) { return 2.0 + 3.0 * x.z; });

  // source_rate * joules_per_energy_unit = 1; the 2.0 cross-sectional area scales volume and power.
  const auto cell_data = coupling::MakeCellPowerData(*field, 0, 2.0, 0.5, 2.0);
  EXPECT_EQ(cell_data.GetDimension(), 1);
  ASSERT_EQ(cell_data.GetCells().size(), grid->GetLocalCellCount());
  for (const auto& cell : cell_data.GetCells())
  {
    EXPECT_EQ(cell.block_id, 7);
    EXPECT_NEAR(cell.volume, 2.0, 1.0e-13);
    EXPECT_NEAR(cell.power, 2.0 * (2.0 + 3.0 * cell.centroid.z), 1.0e-13);
  }
  // 2 * integral_0^4 (2 + 3z) dz = 2 * 32.
  EXPECT_NEAR(GlobalPower(cell_data), 64.0, 1.0e-12);
}

TEST(PreciceCellPower, IntegratesOwnedLinearField2D)
{
  auto grid = BuildSquareMesh(2.0, 2, 0.0);
  const auto field = MakeNodalField(grid, [](const Vector3& x) { return 1.0 + x.x + 2.0 * x.y; });

  const auto cell_data = coupling::MakeCellPowerData(*field, 0, 1.0, 1.0, 0.5);
  EXPECT_EQ(cell_data.GetDimension(), 2);
  for (const auto& cell : cell_data.GetCells())
  {
    EXPECT_NEAR(cell.volume, 0.5, 1.0e-13);
    EXPECT_NEAR(cell.power, 0.5 * (1.0 + cell.centroid.x + 2.0 * cell.centroid.y), 1.0e-13);
  }
  // 0.5 * integral over [0,2]^2 of (1 + x + 2y) = 0.5 * (4 + 4 + 8).
  EXPECT_NEAR(GlobalPower(cell_data), 8.0, 1.0e-12);
}

TEST(PreciceCellPower, IntegratesOwnedTrilinearField3D)
{
  auto grid = BuildBoxMesh(2.0, 2, 0.0);
  const auto field = MakeNodalField(grid, [](const Vector3& x) { return x.x * x.y * x.z; });

  const auto cell_data = coupling::MakeCellPowerData(*field, 0, 1.0, 1.0);
  EXPECT_EQ(cell_data.GetDimension(), 3);
  for (const auto& cell : cell_data.GetCells())
  {
    EXPECT_NEAR(cell.volume, 1.0, 1.0e-13);
    // The mean of a product of independent linear factors over a box is its centroid value.
    EXPECT_NEAR(cell.power, cell.centroid.x * cell.centroid.y * cell.centroid.z, 1.0e-13);
  }
  EXPECT_NEAR(GlobalPower(cell_data), 8.0, 1.0e-12);
}

TEST(PreciceCellPower, SelectsBlocks)
{
  auto grid = BuildLineMesh(4.0, 4, 0.0);
  for (const auto& cell : grid->GetLocalCells())
    cell->block_id = cell->centroid.z < 2.0 ? 1 : 2;
  const auto field = MakeNodalField(grid, [](const Vector3& x) { return 2.0 + 3.0 * x.z; });

  const auto cell_data = coupling::MakeCellPowerData(*field, 0, 1.0, 1.0, 1.0, {2, 2});
  for (const auto& cell : cell_data.GetCells())
    EXPECT_EQ(cell.block_id, 2);
  // integral_2^4 (2 + 3z) dz.
  EXPECT_NEAR(GlobalPower(cell_data), 22.0, 1.0e-12);

  EXPECT_THROW(coupling::MakeCellPowerData(*field, 0, 1.0, 1.0, 1.0, {2, 3}),
               std::invalid_argument);
}

TEST(PreciceCellPower, RejectsInvalidSettings)
{
  auto line = BuildLineMesh(1.0, static_cast<unsigned int>(mpi_comm.size()), 0.0);
  const auto line_field = MakeNodalField(line, [](const Vector3&) { return 1.0; });

  EXPECT_THROW(coupling::MakeCellPowerData(*line_field, 1, 1.0, 1.0, 1.0), std::out_of_range);
  EXPECT_THROW(coupling::MakeCellPowerData(*line_field, 0, -1.0, 1.0, 1.0), std::invalid_argument);
  EXPECT_THROW(coupling::MakeCellPowerData(*line_field, 0, 1.0, 0.0, 1.0), std::invalid_argument);
  EXPECT_THROW(
    coupling::MakeCellPowerData(*line_field, 0, std::numeric_limits<double>::quiet_NaN(), 1.0, 1.0),
    std::invalid_argument);
  EXPECT_THROW(coupling::MakeCellPowerData(*line_field, 0, 1.0, 1.0), std::invalid_argument);
  EXPECT_THROW(coupling::MakeCellPowerData(*line_field, 0, 1.0, 1.0, 0.0), std::invalid_argument);
  EXPECT_THROW(
    coupling::MakeCellPowerData(*line_field, 0, 1.0, 1.0, std::numeric_limits<double>::infinity()),
    std::invalid_argument);

  auto box = BuildBoxMesh(1.0, 2, 0.0);
  const auto box_field = MakeNodalField(box, [](const Vector3&) { return 1.0; });
  EXPECT_THROW(coupling::MakeCellPowerData(*box_field, 0, 1.0, 1.0, 1.0), std::invalid_argument);

  auto rz = BuildSquareMesh(1.0, 2, 0.0);
  rz->SetCoordinateSystem(CoordinateSystemType::CYLINDRICAL);
  const auto rz_field = MakeNodalField(rz, [](const Vector3&) { return 1.0; });
  EXPECT_THROW(coupling::MakeCellPowerData(*rz_field, 0, 1.0, 1.0, 1.0), std::invalid_argument);
}

TEST(PreciceCellPower, PropagatesRankLocalFailures)
{
  auto grid = BuildLineMesh(4.0, 4, 0.0);
  const auto field = MakeNodalField(grid, [](const Vector3&) { return 1.0; });
  // Corrupt one value on rank 0 only; every rank must throw rather than continue alone.
  if (mpi_comm.rank() == 0)
    field->GetLocalFieldVector()[0] = std::numeric_limits<double>::quiet_NaN();
  EXPECT_THROW(coupling::MakeCellPowerData(*field, 0, 1.0, 1.0, 1.0), std::runtime_error);
}

TEST(PreciceCellPower, ValidatesCellRecords)
{
  const coupling::OwnedCellPower valid{3, 1, Vector3(0.0, 0.0, 0.0), 2.0, 4.0};
  EXPECT_DOUBLE_EQ(coupling::CellPowerData(3, {valid}).GetLocalPower(), 4.0);
  EXPECT_THROW(coupling::CellPowerData(3, {valid, valid}), std::invalid_argument);
  EXPECT_THROW(coupling::CellPowerData(0, {valid}), std::invalid_argument);
  EXPECT_THROW(coupling::CellPowerData(4, {valid}), std::invalid_argument);

  auto invalid_volume = valid;
  invalid_volume.volume = 0.0;
  EXPECT_THROW(coupling::CellPowerData(3, {invalid_volume}), std::invalid_argument);

  auto invalid_power = valid;
  invalid_power.power = std::numeric_limits<double>::infinity();
  EXPECT_THROW(coupling::CellPowerData(3, {invalid_power}), std::invalid_argument);
}

TEST(PreciceCellPower, ValidatesTransferPowerCollectively)
{
  const auto rank_id = static_cast<std::uint64_t>(mpi_comm.rank());
  coupling::CellPowerData valid(3, {{rank_id, 1, Vector3(), 1.0, 2.0}});
  EXPECT_DOUBLE_EQ(coupling::GetValidatedGlobalPower(valid), 2.0 * mpi_comm.size());

  coupling::CellPowerData negative(3, {{rank_id, 1, Vector3(), 1.0, -1.0}});
  EXPECT_THROW(coupling::GetValidatedGlobalPower(negative), std::invalid_argument);

  coupling::CellPowerData empty(3, {});
  EXPECT_THROW(coupling::GetValidatedGlobalPower(empty), std::invalid_argument);

  if (mpi_comm.size() > 1)
  {
    coupling::CellPowerData mixed(mpi_comm.rank() == 0 ? 2 : 3,
                                  {{rank_id, 1, Vector3(), 1.0, 2.0}});
    EXPECT_THROW(coupling::GetValidatedGlobalPower(mixed), std::invalid_argument);
  }
}

TEST(PreciceCellPower, ConservativelyAveragesGroups)
{
  // Every rank holds the same cells, so each group's global density equals its local density.
  coupling::CellPowerData data(3,
                               {{0, 1, Vector3(0.0, 0.0, 0.0), 1.0, -1.0},
                                {1, 1, Vector3(1.0, 0.0, 0.0), 3.0, 5.0},
                                {2, 2, Vector3(2.0, 0.0, 0.0), 2.0, 6.0}});
  const double original_global_power = GlobalPower(data);
  std::vector<std::int64_t> groups(data.GetCells().size());
  for (std::size_t i = 0; i < groups.size(); ++i)
    groups[i] = data.GetCells()[i].block_id;

  data.ConservativelyAverageByGroup(groups);
  EXPECT_NEAR(GlobalPower(data), original_global_power, 1.0e-13 * mpi_comm.size());
  // Group 1 has power 4 in volume 4 and group 2 has power 6 in volume 2.
  const auto& cells = data.GetCells();
  EXPECT_NEAR(cells[0].power, 1.0, 1.0e-14);
  EXPECT_NEAR(cells[1].power, 3.0, 1.0e-14);
  EXPECT_NEAR(cells[2].power, 6.0, 1.0e-14);

  coupling::CellPowerData negative_data(3, {{0, 1, Vector3(0.0, 0.0, 0.0), 1.0, -1.0}});
  EXPECT_THROW(negative_data.ConservativelyAverageByGroup({1}), std::runtime_error);
  EXPECT_THROW(negative_data.ConservativelyAverageByGroup({}), std::invalid_argument);
}

TEST(PreciceCellPower, AveragesGroupsAcrossRanks)
{
  // Rank 0 owns no cells. On the others, rank r owns one cell of volume r and power r * r, all in
  // group 5, so the global density is sum(r^2) / sum(r).
  std::vector<coupling::OwnedCellPower> cells;
  std::vector<std::int64_t> groups;
  const int rank = mpi_comm.rank();
  if (rank > 0)
  {
    const double r = rank;
    cells.push_back({static_cast<std::uint64_t>(rank), 1, Vector3(r, 0.0, 0.0), r, r * r});
    groups.push_back(5);
  }
  coupling::CellPowerData data(3, std::move(cells));
  const double original_global_power = GlobalPower(data);
  data.ConservativelyAverageByGroup(groups);
  EXPECT_NEAR(GlobalPower(data), original_global_power, 1.0e-12);

  double sum_r = 0.0;
  double sum_r2 = 0.0;
  for (int r = 1; r < mpi_comm.size(); ++r)
  {
    sum_r += r;
    sum_r2 += static_cast<double>(r) * r;
  }
  if (rank > 0)
    EXPECT_NEAR(data.GetCells()[0].power, rank * sum_r2 / sum_r, 1.0e-12);
}

#endif // OPENSN_WITH_PRECICE
