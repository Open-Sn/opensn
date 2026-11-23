// SPDX-FileCopyrightText: 2026 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "test/unit/common/mesh_builders.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/discrete_ordinates_problem.h"
#include "framework/math/linear_solver/linear_solver.h"
#include "framework/materials/multi_group_xs/multi_group_xs.h"
#include "framework/math/quadratures/angular/product_quadrature.h"
#include "framework/parameters/parameter_block.h"
#include "framework/runtime.h"
#include <gtest/gtest.h>
#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstddef>
#include <filesystem>
#include <fstream>
#include <memory>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

using namespace opensn;

namespace
{

/// Three groups with within-group, down-, and upscatter, so every groupset solve builds scattering
/// sources from the other groupset as well as its own.
std::shared_ptr<MultiGroupXS>
MakeScatteringXS()
{
  const auto now = std::chrono::steady_clock::now().time_since_epoch().count();
  const auto path = std::filesystem::temp_directory_path() /
                    ("opensn_groupset_solve_q_r" + std::to_string(mpi_comm.rank()) + "_" +
                     std::to_string(now) + ".xs");
  {
    std::ofstream out(path);
    out << "NUM_GROUPS 3\nNUM_MOMENTS 2\n"
           "SIGMA_T_BEGIN\n0 1.0\n1 1.5\n2 2.0\nSIGMA_T_END\n"
           "TRANSFER_MOMENTS_BEGIN\n"
           "M_GFROM_GTO_VAL 0 0 0 0.5\n"
           "M_GFROM_GTO_VAL 0 0 1 0.2\n"
           "M_GFROM_GTO_VAL 0 1 1 0.7\n"
           "M_GFROM_GTO_VAL 0 1 2 0.3\n"
           "M_GFROM_GTO_VAL 0 2 1 0.1\n"
           "M_GFROM_GTO_VAL 0 2 2 1.2\n"
           "M_GFROM_GTO_VAL 1 0 0 0.1\n"
           "M_GFROM_GTO_VAL 1 1 1 0.1\n"
           "M_GFROM_GTO_VAL 1 2 2 0.2\n"
           "TRANSFER_MOMENTS_END\n";
  }
  auto xs = std::make_shared<MultiGroupXS>(MultiGroupXS::LoadFromOpenSn(path.string()));
  std::filesystem::remove(path);
  return xs;
}

/// A 2D P1 problem with groupsets (0) and (1, 2) solved by `method`. With `reflecting`, the
/// opposing reflecting boundaries xmin and xmax give the groupset solves delayed angular unknowns;
/// a single reflecting boundary would not. The first groupset's solver then owns its vectors, and
/// the second's, whose block is last, views the flux storage and its tail.
std::shared_ptr<DiscreteOrdinatesProblem>
MakeProblem(const std::string& method, int max_iterations, bool reflecting)
{
  auto grid = BuildSquareMesh(1.0, 4, 0.0);
  grid->SetUniformBlockID(0);

  const auto quadrature =
    std::shared_ptr<AngularQuadrature>(std::make_shared<GLCProductQuadrature2DXY>(2, 4, 1));
  ParameterBlock groupsets("groupsets");
  const std::vector<std::vector<int>> ranges = {{0, 0}, {1, 2}};
  for (std::size_t gs = 0; gs < ranges.size(); ++gs)
  {
    ParameterBlock groupset(std::to_string(gs));
    groupset.AddParameter(ParameterBlock("groups_from_to", ranges[gs]));
    groupset.AddParameter("angular_quadrature", quadrature);
    groupset.AddParameter("inner_linear_method", method);
    groupset.AddParameter("l_max_its", max_iterations);
    groupset.AddParameter("l_abs_tol", 1.0e-10);
    groupsets.AddParameter(groupset);
  }
  groupsets.ChangeToArray();

  ParameterBlock xs_entry("0");
  xs_entry.AddParameter(ParameterBlock("block_ids", std::vector<int>{0}));
  xs_entry.AddParameter("xs", MakeScatteringXS());
  ParameterBlock xs_map("xs_map");
  xs_map.AddParameter(xs_entry);
  xs_map.ChangeToArray();

  ParameterBlock params;
  params.AddParameter("mesh", grid);
  params.AddParameter("num_groups", 3);
  params.AddParameter(groupsets);
  params.AddParameter(xs_map);
  if (reflecting)
  {
    ParameterBlock boundaries("boundary_conditions");
    for (const std::string name : {"xmin", "xmax"})
    {
      ParameterBlock boundary(std::to_string(boundaries.GetNumParameters()));
      boundary.AddParameter("name", name);
      boundary.AddParameter("type", std::string("reflecting"));
      boundaries.AddParameter(boundary);
    }
    boundaries.ChangeToArray();
    params.AddParameter(boundaries);
  }
  return DiscreteOrdinatesProblem::Create(params);
}

/// Inner method, maximum iterations (1 makes PETSc Richardson a single sweep), reflecting boundary.
using SolveCase = std::tuple<std::string, int, bool>;

class GroupsetSolveSourceMomentsTest : public testing::TestWithParam<SolveCase>
{
};

} // namespace

// The within-group solvers build their sources in scratch, and the across-groupset solve relies on
// them leaving the source moments unchanged, including the sweep that rebuilds the angular flux
// after a Krylov solve.
TEST_P(GroupsetSolveSourceMomentsTest, LeavesSourceMomentsUnchanged)
{
  const auto& [method, max_iterations, reflecting] = GetParam();
  auto problem = MakeProblem(method, max_iterations, reflecting);

  auto& q_moments = problem->GetQMomentsLocal();
  std::vector<double> q_before(q_moments.Size());
  for (std::size_t k = 0; k < q_before.size(); ++k)
    q_before[k] = 1.0 + 0.01 * static_cast<double>(k % 97);
  q_moments.Assign(q_before);

  for (std::size_t gs = 0; gs < problem->GetNumWGSSolvers(); ++gs)
  {
    auto solver = problem->GetWGSSolver(gs);
    solver->Setup();
    solver->Solve();

    const auto q_after = std::as_const(q_moments).Data();
    ASSERT_EQ(q_after.size(), q_before.size());
    for (std::size_t k = 0; k < q_before.size(); ++k)
      ASSERT_EQ(q_after[k], q_before[k]) << "groupset " << gs << ", entry " << k;
  }

  // The solves did work: every groupset has a positive, finite flux.
  const auto& groupsets = problem->GetGroupsets();
  for (std::size_t gs = 0; gs < groupsets.size(); ++gs)
  {
    const auto block =
      std::as_const(problem->GetPhiNewLocal()).Block(problem->GetFluxBlock(groupsets[gs]));
    double local_max = 0.0;
    for (const double value : block)
    {
      ASSERT_TRUE(std::isfinite(value));
      local_max = std::max(local_max, value);
    }
    double global_max = 0.0;
    mpi_comm.all_reduce(local_max, global_max, mpi::op::max<double>());
    EXPECT_GT(global_max, 0.0) << "groupset " << gs;
  }
}

INSTANTIATE_TEST_SUITE_P(InnerMethods,
                         GroupsetSolveSourceMomentsTest,
                         testing::Values(SolveCase{"petsc_gmres", 100, false},
                                         SolveCase{"petsc_gmres", 100, true},
                                         SolveCase{"petsc_bicgstab", 100, true},
                                         SolveCase{"petsc_richardson", 1, false},
                                         SolveCase{"petsc_richardson", 1, true},
                                         SolveCase{"petsc_richardson", 100, true},
                                         SolveCase{"classic_richardson", 100, false},
                                         SolveCase{"classic_richardson", 100, true}),
                         [](const testing::TestParamInfo<SolveCase>& info)
                         {
                           // Not a structured binding: its commas would split the macro arguments.
                           return std::get<0>(info.param) + "_" +
                                  std::to_string(std::get<1>(info.param)) +
                                  (std::get<2>(info.param) ? "_reflecting" : "_vacuum");
                         });
