// SPDX-FileCopyrightText: 2026 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "framework/field_functions/field_function_grid_based.h"
#include "framework/field_functions/interpolation/ffinter_volume.h"
#include "framework/math/spatial_discretization/finite_element/piecewise_linear/piecewise_linear_discontinuous.h"
#include "framework/mesh/logical_volume/rpp_logical_volume.h"
#include "test/unit/common/mesh_builders.h"
#include "gtest/gtest.h"

using namespace opensn;

TEST(FieldFunctionInterpolationVolume, TransformedExtrema)
{
  auto grid = BuildLineMesh(4.0, 4, 0.0);
  std::shared_ptr<SpatialDiscretization> sdm = PieceWiseLinearDiscontinuous::New(grid);
  auto ff = std::make_shared<FieldFunctionGridBased>("linear", sdm, Unknown(UnknownType::SCALAR));
  const auto& unknowns = ff->GetUnknownManager();
  std::vector<double> values(sdm->GetNumLocalDOFs(unknowns));
  for (auto& cell : grid->local_cells)
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
