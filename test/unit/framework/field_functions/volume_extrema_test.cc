// SPDX-FileCopyrightText: 2026 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "framework/field_functions/field_function_grid_based.h"
#include "framework/field_functions/interpolation/ffinter_volume.h"
#include "framework/math/spatial_discretization/finite_element/piecewise_linear/piecewise_linear_discontinuous.h"
#include "framework/mesh/logical_volume/rpp_logical_volume.h"
#include "test/unit/common/mesh_builders.h"
#include "gtest/gtest.h"
#include <cmath>
#include <filesystem>

using namespace opensn;

TEST(FieldFunctionInterpolationVolume, TransformedExtrema)
{
  auto grid = BuildLineMesh(4.0, 4, 0.0);
  std::shared_ptr<SpatialDiscretization> sdm = PieceWiseLinearDiscontinuous::New(grid);
  auto ff = std::make_shared<FieldFunctionGridBased>("linear", sdm, Unknown(UnknownType::SCALAR));
  const auto& unknowns = ff->GetUnknownManager();
  std::vector<double> values(sdm->GetNumLocalDOFs(unknowns));
  for (auto& cell : grid->GetLocalCells())
  {
    cell.block_id = 3;
    const auto& nodes = sdm->GetCellNodeLocations(cell);
    for (size_t i = 0; i < nodes.size(); ++i)
      values[sdm->MapDOFLocal(cell, i, unknowns, 0, 0)] = nodes[i].z;
  }
  ff->UpdateFieldVector(values);

  auto params = RPPLogicalVolume::GetInputParameters();
  ParameterBlock volume_params;
  volume_params.AddParameter("infx", true);
  volume_params.AddParameter("infy", true);
  volume_params.AddParameter("infz", true);
  params.AssignParameters(volume_params);
  FieldFunctionInterpolationVolume interp;
  interp.SetLogicalVolume(std::make_shared<RPPLogicalVolume>(params));
  interp.AddFieldFunction(ff);
  // Reverses the extrema and shifts all transformed values above the raw range.
  interp.SetOperationFunction([](double value, unsigned int block_id)
                              { return 10.0 + block_id - 2.0 * value; });

  for (const auto& [operation, expected] :
       {std::pair{FieldFunctionInterpolationOperation::OP_MIN_FUNC, 5.0},
        std::pair{FieldFunctionInterpolationOperation::OP_MAX_FUNC, 13.0},
        std::pair{FieldFunctionInterpolationOperation::OP_MIN, 0.0},
        std::pair{FieldFunctionInterpolationOperation::OP_MAX, 4.0}})
  {
    interp.SetOperationType(operation);
    interp.Execute();
    EXPECT_NEAR(interp.GetValue(), expected, 1.0e-12);
  }

  // Restrict to one cell so an MPI run also exercises ranks with no selected cells.
  auto restricted_params = RPPLogicalVolume::GetInputParameters();
  ParameterBlock restricted_volume;
  restricted_volume.AddParameter("infx", true);
  restricted_volume.AddParameter("infy", true);
  restricted_params.AssignParameters(restricted_volume);
  interp.SetLogicalVolume(std::make_shared<RPPLogicalVolume>(restricted_params));
  interp.SetOperationType(FieldFunctionInterpolationOperation::OP_MIN_FUNC);
  interp.Execute();
  EXPECT_NEAR(interp.GetValue(), 11.0, 1.0e-12);
  interp.SetOperationType(FieldFunctionInterpolationOperation::OP_MAX_FUNC);
  interp.Execute();
  EXPECT_NEAR(interp.GetValue(), 13.0, 1.0e-12);
}

TEST(FieldFunctionGridBased, BatchedPointValues)
{
  auto grid = BuildLineMesh(4.0, 4, 0.0);
  std::shared_ptr<SpatialDiscretization> sdm = PieceWiseLinearDiscontinuous::New(grid);
  auto linear =
    std::make_shared<FieldFunctionGridBased>("linear", sdm, Unknown(UnknownType::SCALAR));
  auto shifted =
    std::make_shared<FieldFunctionGridBased>("shifted", sdm, Unknown(UnknownType::SCALAR));
  const auto& unknowns = linear->GetUnknownManager();
  std::vector<double> linear_values(sdm->GetNumLocalDOFs(unknowns));
  std::vector<double> shifted_values(linear_values.size());
  for (const auto& cell : grid->GetLocalCells())
  {
    const auto& nodes = sdm->GetCellNodeLocations(cell);
    for (size_t i = 0; i < nodes.size(); ++i)
    {
      const auto dof = sdm->MapDOFLocal(cell, i, unknowns, 0, 0);
      linear_values[dof] = nodes[i].z;
      shifted_values[dof] = 2.0 * nodes[i].z + 3.0;
    }
  }
  linear->UpdateFieldVector(linear_values);
  shifted->UpdateFieldVector(shifted_values);

  const std::vector<Vector3> points{{0.0, 0.0, 0.25}, {0.0, 0.0, 1.0}, {0.0, 0.0, 2.75}};
  const FieldFunctionGridBased::FFList fields{linear, shifted};
  const auto values = FieldFunctionGridBased::GetPointValues(fields, points);

  ASSERT_EQ(values.size(), 2);
  ASSERT_EQ(values[0].size(), points.size());
  for (size_t p = 0; p < points.size(); ++p)
  {
    ASSERT_EQ(values[0][p].size(), 1);
    EXPECT_NEAR(values[0][p][0], points[p].z, 1.0e-12);
    EXPECT_NEAR(values[1][p][0], 2.0 * points[p].z + 3.0, 1.0e-12);
    EXPECT_NEAR(values[0][p][0], linear->GetPointValue(points[p])[0], 1.0e-12);
    EXPECT_NEAR(values[1][p][0], shifted->GetPointValue(points[p])[0], 1.0e-12);
  }
}

// The binned point location of GetPointValues must find the same cells as the brute-force search
// of GetPointValue. On the vertices, edges, and faces of the cells, a point lies in several cells
// and gets the average of their values. A linear field is exact on every cell, so both searches
// must also return the field itself. A field that is constant on each cell and differs from cell
// to cell changes the average at such a point when a search misses one of its cells.
TEST(FieldFunctionGridBased, BatchedPointValues3D)
{
  const auto tet_file = std::filesystem::path(OPENSN_TEST_ROOT) / "assets" / "mesh" / "cube3.2.msh";
  // A domain 1e-6 thick along z is thicker than the bin padding, so the bins drop that axis and
  // spread over x and y.
  std::vector<double> thin_nodes(7);
  for (size_t i = 0; i < thin_nodes.size(); ++i)
    thin_nodes[i] = static_cast<double>(i) / 6.0;
  const std::vector<std::pair<std::shared_ptr<Mesh>, Vector3>> meshes = {
    {BuildMeshFromFile(tet_file), Vector3(0.032, 0.032, 0.032)},
    {BuildBoxMesh(1.0, 4, 0.0), Vector3(1.0, 1.0, 1.0)},
    {BuildOrthogonalMesh({thin_nodes, thin_nodes, {0.0, 1.0e-6}}), Vector3(1.0, 1.0, 1.0e-6)}};
  for (const auto& [grid, size] : meshes)
  {
    std::shared_ptr<SpatialDiscretization> sdm = PieceWiseLinearDiscontinuous::New(grid);
    auto ff = std::make_shared<FieldFunctionGridBased>("linear", sdm, Unknown(UnknownType::SCALAR));
    const auto field = [size](const Vector3& p)
    { return 1.0 + p.x / size.x + 2.0 * p.y / size.y + 3.0 * p.z / size.z; };
    const auto& unknowns = ff->GetUnknownManager();
    std::vector<double> values(sdm->GetNumLocalDOFs(unknowns));
    for (const auto& cell : grid->GetLocalCells())
    {
      const auto& nodes = sdm->GetCellNodeLocations(cell);
      for (size_t i = 0; i < nodes.size(); ++i)
        values[sdm->MapDOFLocal(cell, i, unknowns, 0, 0)] = field(nodes[i]);
    }
    ff->UpdateFieldVector(values);

    auto jump = std::make_shared<FieldFunctionGridBased>("jump", sdm, Unknown(UnknownType::SCALAR));
    std::vector<double> jump_values(sdm->GetNumLocalDOFs(unknowns));
    for (const auto& cell : grid->GetLocalCells())
      for (size_t i = 0; i < sdm->GetCellNumNodes(cell); ++i)
        jump_values[sdm->MapDOFLocal(cell, i, unknowns, 0, 0)] =
          1.0 + static_cast<double>(cell.global_id);
    jump->UpdateFieldVector(jump_values);

    // A lattice that includes the domain boundary and, on the box mesh, the vertices, edges, and
    // faces of the cells, and generic points. Every rank uses the same points.
    std::vector<Vector3> points;
    const int n = 8;
    for (int k = 0; k <= n; ++k)
      for (int j = 0; j <= n; ++j)
        for (int i = 0; i <= n; ++i)
          points.emplace_back(size.x * i / n, size.y * j / n, size.z * k / n);
    // Points of the R3 low-discrepancy sequence, spread over the domain.
    const double g = 1.22074408460575947536;
    for (int p = 1; p <= 200; ++p)
    {
      const double x = std::fmod(0.5 + p / g, 1.0);
      const double y = std::fmod(0.5 + p / (g * g), 1.0);
      const double z = std::fmod(0.5 + p / (g * g * g), 1.0);
      points.emplace_back(size.x * x, size.y * y, size.z * z);
    }

    const auto batched = FieldFunctionGridBased::GetPointValues({ff, jump}, points);
    ASSERT_EQ(batched.size(), 2);
    ASSERT_EQ(batched[0].size(), points.size());
    ASSERT_EQ(batched[1].size(), points.size());
    for (size_t p = 0; p < points.size(); ++p)
    {
      const double expected = field(points[p]);
      EXPECT_NEAR(batched[0][p][0], expected, 1.0e-12 * expected);
      EXPECT_NEAR(batched[0][p][0], ff->GetPointValue(points[p])[0], 1.0e-14 * expected);
      const double brute_force_jump = jump->GetPointValue(points[p])[0];
      EXPECT_NEAR(batched[1][p][0], brute_force_jump, 1.0e-14 * brute_force_jump);
    }
  }
}
