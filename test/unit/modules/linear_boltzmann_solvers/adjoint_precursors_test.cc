// SPDX-FileCopyrightText: 2026 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "test/unit/common/mesh_builders.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/discrete_ordinates_problem.h"
#include "modules/linear_boltzmann_solvers/lbs_problem/compute/lbs_compute.h"
#include "framework/materials/multi_group_xs/multi_group_xs.h"
#include "framework/math/quadratures/angular/product_quadrature.h"
#include "framework/parameters/parameter_block.h"
#include <gtest/gtest.h>
#include <memory>
#include <string>
#include <vector>

using namespace opensn;

namespace
{

/// A one-group 2D problem with one delayed-neutron precursor family and use_precursors=true.
std::shared_ptr<DiscreteOrdinatesProblem>
MakeDelayedProblem(bool adjoint)
{
  auto grid = BuildSquareMesh(1.0, 4, 0.0);
  grid->SetUniformBlockID(0);

  ParameterBlock groupset("0");
  groupset.AddParameter(ParameterBlock("groups_from_to", std::vector<int>{0, 0}));
  groupset.AddParameter(
    "angular_quadrature",
    std::shared_ptr<AngularQuadrature>(std::make_shared<GLCProductQuadrature2DXY>(2, 4, 0)));
  ParameterBlock groupsets("groupsets");
  groupsets.AddParameter(groupset);
  groupsets.ChangeToArray();

  auto xs = std::make_shared<MultiGroupXS>(MultiGroupXS::LoadFromOpenSn(
    std::string(OPENSN_TEST_ROOT) + "/assets/xs/xs1g_delayed_crit_1p.cxs"));
  ParameterBlock xs_entry("0");
  xs_entry.AddParameter(ParameterBlock("block_ids", std::vector<int>{0}));
  xs_entry.AddParameter("xs", xs);
  ParameterBlock xs_map("xs_map");
  xs_map.AddParameter(xs_entry);
  xs_map.ChangeToArray();

  ParameterBlock options("options");
  options.AddParameter("adjoint", adjoint);
  options.AddParameter("use_precursors", true);

  ParameterBlock params;
  params.AddParameter("mesh", grid);
  params.AddParameter("num_groups", 1);
  params.AddParameter(groupsets);
  params.AddParameter(xs_map);
  params.AddParameter(options);
  return DiscreteOrdinatesProblem::Create(params);
}

} // namespace

// Equilibrium precursor concentrations follow from a forward flux, C = gamma / lambda * nu_d
// sigma_f phi. An adjoint flux has none, so ComputePrecursors leaves them zero in adjoint mode.
TEST(AdjointPrecursorsTest, ComputedOnlyFromForwardFlux)
{
  // xs1g_delayed_crit_1p.cxs: lambda = 0.08, nu_d = 0.013, sigma_f = 0.15. The yield of its single
  // family is normalized to gamma = 1.
  const double expected = 1.0 / 0.08 * 0.013 * 0.15;
  for (const bool adjoint : {false, true})
  {
    auto problem = MakeDelayedProblem(adjoint);
    problem->SetPhiNewFrom(std::vector<double>(problem->GetPhiNewLocal().size(), 1.0));
    auto& precursors = problem->GetPrecursorsNewLocal();
    ASSERT_FALSE(precursors.empty());
    precursors.assign(precursors.size(), -1.0);

    ComputePrecursors(*problem);
    for (const auto value : problem->GetPrecursorsNewLocal())
    {
      if (adjoint)
        EXPECT_EQ(value, 0.0);
      else
        EXPECT_DOUBLE_EQ(value, expected);
    }
  }
}
