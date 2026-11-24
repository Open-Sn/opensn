// SPDX-FileCopyrightText: 2026 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "test/unit/common/mesh_builders.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/discrete_ordinates_problem.h"
#include "modules/linear_boltzmann_solvers/lbs_problem/compute/lbs_compute.h"
#include "framework/materials/multi_group_xs/multi_group_xs.h"
#include "framework/math/quadratures/angular/product_quadrature.h"
#include "framework/math/spatial_discretization/spatial_discretization.h"
#include "framework/mesh/logical_volume/rpp_logical_volume.h"
#include "framework/parameters/parameter_block.h"
#include "framework/runtime.h"
#include <gtest/gtest.h>
#include <algorithm>
#include <cmath>
#include <memory>
#include <string>
#include <utility>
#include <vector>

using namespace opensn;

namespace
{

std::pair<std::shared_ptr<DiscreteOrdinatesProblem>, std::shared_ptr<MultiGroupXS>>
MakeProductionProblem(const bool adjoint, const bool delayed)
{
  auto grid = BuildSquareMesh(1.0, 4, 0.0);
  grid->SetUniformBlockID(0);
  ParameterBlock right_half;
  right_half.AddParameter("xmin", 0.5);
  right_half.AddParameter("xmax", 1.0);
  right_half.AddParameter("infy", true);
  right_half.AddParameter("infz", true);
  grid->SetBlockIDFromLogicalVolume(*RPPLogicalVolume::Create(right_half), 1, true);

  ParameterBlock groupset("0");
  groupset.AddParameter(ParameterBlock("groups_from_to", std::vector<int>{0, 1}));
  groupset.AddParameter(
    "angular_quadrature",
    std::shared_ptr<AngularQuadrature>(std::make_shared<GLCProductQuadrature2DXY>(2, 4, 0)));
  ParameterBlock groupsets("groupsets");
  groupsets.AddParameter(groupset);
  groupsets.ChangeToArray();

  ParameterBlock xs_map("xs_map");
  std::shared_ptr<MultiGroupXS> right_xs;
  for (int block = 0; block < 2; ++block)
  {
    const std::string filename =
      block == 0 ? "xs2g_delayed_crit_1p.cxs" : "xs2g_delayed_density_up_1p.cxs";
    auto xs = std::make_shared<MultiGroupXS>(
      MultiGroupXS::LoadFromOpenSn(std::string(OPENSN_TEST_ROOT) + "/assets/xs/" + filename));
    if (block == 1)
      right_xs = xs;
    ParameterBlock entry(std::to_string(block));
    entry.AddParameter(ParameterBlock("block_ids", std::vector<int>{block}));
    entry.AddParameter("xs", xs);
    xs_map.AddParameter(entry);
  }
  xs_map.ChangeToArray();
  ParameterBlock options("options");
  options.AddParameter("adjoint", adjoint);
  options.AddParameter("use_precursors", delayed);
  ParameterBlock params;
  params.AddParameter("mesh", grid);
  params.AddParameter("num_groups", 2);
  params.AddParameter(groupsets);
  params.AddParameter(xs_map);
  params.AddParameter(options);
  return {DiscreteOrdinatesProblem::Create(params), right_xs};
}

} // namespace

TEST(FissionProductionTest, AnalyticAffineFluxWithDelayedAdjointAndScaledXS)
{
  for (const bool adjoint : {false, true})
    for (const bool delayed : {false, true})
    {
      SCOPED_TRACE(::testing::Message() << "adjoint=" << adjoint << ", delayed=" << delayed);
      auto [problem, right_xs] = MakeProductionProblem(adjoint, delayed);
      const auto& sdm = problem->GetSpatialDiscretization();
      const auto& unknowns = problem->GetUnknownManager();
      auto phi = problem->GetPhiNewLocal();
      for (const auto& cell : problem->GetGrid()->GetLocalCells())
      {
        const auto& nodes = sdm.GetCellNodeLocations(*cell);
        for (unsigned int i = 0; i < nodes.size(); ++i)
        {
          const auto& p = nodes[i];
          phi[sdm.MapDOFLocal(*cell, i, unknowns, 0, 0)] = 1.0 + p.x + 2.0 * p.y;
          // Signed scalar nodal values must be integrated, not clipped or discarded.
          phi[sdm.MapDOFLocal(*cell, i, unknowns, 0, 1)] = -0.5 + 0.25 * p.x + 0.5 * p.y;
        }
      }

      // Each half has area 1/2. The exact mean fast/thermal fluxes are
      // (2.25, -0.1875) on the left and (2.75, -0.0625) on the right.
      // sigma_f=(.12,.03), chi_prompt=chi_delayed=(1,0), nu_p=1.987 and nu_d=.013.
      // The right material's density is 1.2 times the left before Scale is applied.
      for (const double scale : {1.0, 1.5, 0.75, 1.0})
      {
        SCOPED_TRACE(::testing::Message() << "scale=" << scale);
        right_xs->Scale(scale);
        const double nu = delayed ? 2.0 : 1.987;
        const double left = adjoint ? 0.15 * 2.25 : 0.12 * 2.25 - 0.03 * 0.1875;
        const double right = adjoint ? 0.15 * 2.75 : 0.12 * 2.75 - 0.03 * 0.0625;
        const double expected = 0.5 * nu * (left + 1.2 * scale * right);
        // Affine PWLD integration is exact. This bound covers only ordinary roundoff in local
        // summation, material coefficients, and MPI reductions on 1, 2, and 4 ranks.
        EXPECT_NEAR(ComputeFissionProduction(*problem, phi), expected, 5.0e-13);
      }
    }
}

TEST(FissionProductionTest, NormalizedCellFissionSourcesAndRMSChange)
{
  auto [problem, right_xs] = MakeProductionProblem(false, false);
  const auto& grid = *problem->GetGrid();

  // Unit flux in both groups. Each of the 16 cells has area 1/16. The production per unit flux of
  // the right material is 1.2 times that of the left, so the normalized source of a cell is
  // 1 / 17.6 on the left and 1.2 / 17.6 on the right.
  auto phi = problem->GetPhiNewLocal();
  std::fill(phi.begin(), phi.end(), 1.0);
  const auto previous = ComputeNormalizedCellFissionSources(*problem, phi);
  ASSERT_EQ(previous.size(), grid.GetLocalCellCount());
  double local_sum = 0.0;
  for (const auto& cell : grid.GetLocalCells())
  {
    const double expected = (cell->block_id == 0 ? 1.0 : 1.2) / 17.6;
    EXPECT_NEAR(previous[cell->local_id], expected, 1.0e-15);
    local_sum += previous[cell->local_id];
  }
  double global_sum = 0.0;
  mpi_comm.all_reduce(local_sum, global_sum, mpi::op::sum<double>());
  EXPECT_NEAR(global_sum, 1.0, 1.0e-14);

  // Doubling the flux on the right gives sources 1 / 27.2 and 2.4 / 27.2. The relative changes are
  // 17.6 / 27.2 - 1 on the 8 left cells and 2 * 17.6 / 27.2 - 1 on the 8 right cells. Scaling the
  // whole flux does not change the normalized sources.
  for (const auto& cell : grid.GetLocalCells())
    if (cell->block_id == 1)
    {
      const auto& transport_view = problem->GetCellTransportViews()[cell->local_id];
      for (int i = 0; i < transport_view.GetNumNodes(); ++i)
        for (unsigned int g = 0; g < 2; ++g)
          phi[transport_view.MapDOF(i, 0, g)] = 2.0;
    }
  for (double& value : phi)
    value *= 3.0;
  const auto current = ComputeNormalizedCellFissionSources(*problem, phi);
  const double left_change = 17.6 / 27.2 - 1.0;
  const double right_change = 2.0 * 17.6 / 27.2 - 1.0;
  const double expected_rms =
    std::sqrt(0.5 * (left_change * left_change + right_change * right_change));
  EXPECT_NEAR(ComputeFissionSourceRMSChange(previous, current), expected_rms, 1.0e-14);
  EXPECT_EQ(ComputeFissionSourceRMSChange(current, current), 0.0);
}

TEST(FissionProductionTest, RMSChangeCountsNonpositivePreviousSources)
{
  // Only the first entry of each rank has a positive previous source. Its relative change is 0.5,
  // and all three entries count, so the RMS is sqrt(0.25 / 3) on any number of ranks.
  const std::vector<double> previous = {1.0, 0.0, -1.0};
  const std::vector<double> current = {1.5, 7.0, 3.0};
  EXPECT_NEAR(ComputeFissionSourceRMSChange(previous, current), std::sqrt(0.25 / 3.0), 1.0e-15);

  // A rank without fissionable cells contributes no cells.
  const std::vector<double> empty;
  const bool first_rank_only = mpi_comm.rank() == 0;
  const auto& local_previous = first_rank_only ? previous : empty;
  const auto& local_current = first_rank_only ? current : empty;
  EXPECT_NEAR(
    ComputeFissionSourceRMSChange(local_previous, local_current), std::sqrt(0.25 / 3.0), 1.0e-15);
}
