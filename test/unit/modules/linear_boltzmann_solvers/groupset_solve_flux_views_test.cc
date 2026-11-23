// SPDX-FileCopyrightText: 2026 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "test/unit/common/mesh_builders.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/discrete_ordinates_problem.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/iterative_methods/wgs_context.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep/angle_aggregation/angle_aggregation.h"
#include "modules/linear_boltzmann_solvers/lbs_problem/flux_moment_layout.h"
#include "modules/linear_boltzmann_solvers/lbs_problem/volumetric_source/volumetric_source.h"
#include "framework/math/linear_solver/linear_solver.h"
#include "framework/materials/multi_group_xs/multi_group_xs.h"
#include "framework/math/quadratures/angular/product_quadrature.h"
#include "framework/parameters/parameter_block.h"
#include "framework/runtime.h"
#include <gtest/gtest.h>
#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <memory>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

using namespace opensn;

/*
 * The within-groupset solvers use the problem's scalar flux storage as their solution and
 * right-hand-side vectors: always without delayed angular unknowns, and with them when the
 * groupset's block is the last one, whose delayed unknowns go to the tail after it. Other
 * groupsets with delayed unknowns own their vectors. These tests check the solutions of both kinds
 * against an exact solution.
 *
 * Reference: a homogeneous square with all four boundaries reflecting and a uniform isotropic
 * source is an infinite medium. The exact angular flux is isotropic and the scalar flux is
 * spatially flat, solving sigma_t[g] phi[g] - sum_g' sigma_s0[g' -> g] phi[g'] = q[g]. The
 * spatial and angular discretizations reproduce a flat isotropic solution exactly, so the discrete
 * solution differs from it only by the iteration error. The opposing reflecting boundaries give
 * every groupset delayed angular unknowns.
 */

namespace
{

constexpr unsigned int kNumGroups = 3;
constexpr std::array<double, kNumGroups> kSigmaT = {1.0, 1.5, 2.0};
constexpr std::array<double, kNumGroups> kSource = {1.0, 0.5, 0.25};

struct Transfer
{
  unsigned int from;
  unsigned int to;
  double value;
};

/// Isotropic scattering with down- and upscatter. Each group scatters less than it removes.
constexpr std::array<Transfer, 6> kTransfers = {
  {{0, 0, 0.5}, {0, 1, 0.2}, {1, 1, 0.7}, {1, 2, 0.3}, {2, 1, 0.1}, {2, 2, 1.2}}};

/// The solvers stop at a relative residual of kSolverTolerance. The operator I - DLinvMS has a
/// condition number of a few for these cross sections, so the error stays far below kAcceptance.
constexpr double kSolverTolerance = 1.0e-10;
constexpr double kAcceptance = 1.0e-8;

std::shared_ptr<MultiGroupXS>
MakeXS()
{
  const auto now = std::chrono::steady_clock::now().time_since_epoch().count();
  const auto path = std::filesystem::temp_directory_path() /
                    ("opensn_groupset_flux_views_r" + std::to_string(mpi_comm.rank()) + "_" +
                     std::to_string(now) + ".xs");
  {
    std::ofstream out(path);
    out << "NUM_GROUPS " << kNumGroups << "\nNUM_MOMENTS 1\nSIGMA_T_BEGIN\n";
    for (unsigned int g = 0; g < kNumGroups; ++g)
      out << g << " " << kSigmaT[g] << "\n";
    out << "SIGMA_T_END\nTRANSFER_MOMENTS_BEGIN\n";
    for (const auto& transfer : kTransfers)
      out << "M_GFROM_GTO_VAL 0 " << transfer.from << " " << transfer.to << " " << transfer.value
          << "\n";
    out << "TRANSFER_MOMENTS_END\n";
  }
  auto xs = std::make_shared<MultiGroupXS>(MultiGroupXS::LoadFromOpenSn(path.string()));
  std::filesystem::remove(path);
  return xs;
}

/// Exact infinite-medium scalar flux: the solution of the 3x3 balance system, by Gaussian
/// elimination with partial pivoting.
std::array<double, kNumGroups>
ExactFlux()
{
  std::array<std::array<double, kNumGroups + 1>, kNumGroups> a{};
  for (unsigned int g = 0; g < kNumGroups; ++g)
  {
    a[g][g] = kSigmaT[g];
    a[g][kNumGroups] = kSource[g];
  }
  for (const auto& transfer : kTransfers)
    a[transfer.to][transfer.from] -= transfer.value;
  for (unsigned int k = 0; k < kNumGroups; ++k)
  {
    unsigned int pivot = k;
    for (unsigned int i = k + 1; i < kNumGroups; ++i)
      if (std::abs(a[i][k]) > std::abs(a[pivot][k]))
        pivot = i;
    std::swap(a[k], a[pivot]);
    for (unsigned int i = k + 1; i < kNumGroups; ++i)
    {
      const double factor = a[i][k] / a[k][k];
      for (unsigned int j = k; j <= kNumGroups; ++j)
        a[i][j] -= factor * a[k][j];
    }
  }
  std::array<double, kNumGroups> phi{};
  for (unsigned int k = kNumGroups; k-- > 0;)
  {
    double sum = a[k][kNumGroups];
    for (unsigned int j = k + 1; j < kNumGroups; ++j)
      sum -= a[k][j] * phi[j];
    phi[k] = sum / a[k][k];
  }
  return phi;
}

std::shared_ptr<VolumetricSource>
MakeSource(double scale)
{
  std::vector<double> strength(kSource.begin(), kSource.end());
  for (auto& value : strength)
    value *= scale;
  ParameterBlock params;
  params.AddParameter(ParameterBlock("block_ids", std::vector<int>{0}));
  params.AddParameter(ParameterBlock("group_strength", strength));
  return VolumetricSource::Create(params);
}

/**
 * A 2D P0 problem with a uniform source, solved by `method` with the given groupsets. With
 * `reflecting`, all four boundaries are reflecting, which makes the problem an infinite medium;
 * otherwise they are vacuum, which gives no delayed angular unknowns.
 */
std::shared_ptr<DiscreteOrdinatesProblem>
MakeProblem(const std::vector<std::vector<int>>& ranges, const std::string& method, bool reflecting)
{
  auto grid = BuildSquareMesh(1.0, 4, 0.0);
  grid->SetUniformBlockID(0);

  const auto quadrature =
    std::shared_ptr<AngularQuadrature>(std::make_shared<GLCProductQuadrature2DXY>(2, 4, 0));
  ParameterBlock groupsets("groupsets");
  for (std::size_t gs = 0; gs < ranges.size(); ++gs)
  {
    ParameterBlock groupset(std::to_string(gs));
    groupset.AddParameter(ParameterBlock("groups_from_to", ranges[gs]));
    groupset.AddParameter("angular_quadrature", quadrature);
    groupset.AddParameter("inner_linear_method", method);
    groupset.AddParameter("l_max_its", 1000);
    groupset.AddParameter("l_abs_tol", kSolverTolerance);
    groupsets.AddParameter(groupset);
  }
  groupsets.ChangeToArray();

  ParameterBlock xs_entry("0");
  xs_entry.AddParameter(ParameterBlock("block_ids", std::vector<int>{0}));
  xs_entry.AddParameter("xs", MakeXS());
  ParameterBlock xs_map("xs_map");
  xs_map.AddParameter(xs_entry);
  xs_map.ChangeToArray();

  ParameterBlock sources("volumetric_sources");
  sources.AddParameter("0", MakeSource(1.0));
  sources.ChangeToArray();

  ParameterBlock params;
  params.AddParameter("mesh", grid);
  params.AddParameter("num_groups", static_cast<int>(kNumGroups));
  params.AddParameter(groupsets);
  params.AddParameter(xs_map);
  params.AddParameter(sources);
  if (reflecting)
  {
    ParameterBlock boundaries("boundary_conditions");
    for (const std::string name : {"xmin", "xmax", "ymin", "ymax"})
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

/// Solves the groupsets once, in order. No group scatters into an earlier groupset, so one pass is
/// exact. `release` frees each solver's workspace afterwards, as the across-groupset solver does.
void
SolveGroupsets(DiscreteOrdinatesProblem& problem, bool release)
{
  for (std::size_t gs = 0; gs < problem.GetNumWGSSolvers(); ++gs)
  {
    auto solver = problem.GetWGSSolver(gs);
    solver->Setup();
    solver->Solve();
    if (release)
      solver->ReleaseWorkspace();
  }
}

/// Largest relative error of the scalar flux over local nodes and groups on all ranks.
double
MaxRelativeFluxError(DiscreteOrdinatesProblem& problem)
{
  const auto exact = ExactFlux();
  const auto& layout = problem.GetFluxMomentLayout();
  const auto phi = std::as_const(problem.GetPhiOldLocal()).Data();
  double local_error = 0.0;
  for (std::uint64_t n = 0; n < problem.GetLocalNodeCount(); ++n)
    for (unsigned int g = 0; g < kNumGroups; ++g)
    {
      const double value = phi[layout.MapDOF(n, 0, g)];
      local_error = std::max(
        local_error, std::isfinite(value) ? std::abs(value - exact[g]) / exact[g] : 1.0e300);
    }
  double error = 0.0;
  mpi_comm.all_reduce(local_error, error, mpi::op::max<double>());
  return error;
}

std::uint64_t
GlobalDelayedUnknowns(DiscreteOrdinatesProblem& problem, std::size_t groupset)
{
  return problem.GetGroupsets().at(groupset).angle_agg->GetNumDelayedAngularDOFs().second;
}

using GroupsetRanges = std::vector<std::vector<int>>;
/// Groupset ranges and inner method.
using SolveCase = std::tuple<GroupsetRanges, std::string>;

class GroupsetSolveFluxViewsTest : public testing::TestWithParam<SolveCase>
{
};

} // namespace

TEST(GroupsetSolveFluxViewsReferenceTest, ExactFluxSolvesTheBalance)
{
  // Independent check of the reference, by hand: group 0 has no inscatter, so
  // phi0 = q0 / (sigma_t0 - sigma_s0[0->0]) = 1 / 0.5 = 2. Groups 1 and 2 then solve
  // 0.8 phi1 - 0.1 phi2 = 0.5 + 0.2 * 2 and -0.3 phi1 + 0.8 phi2 = 0.25, which give
  // phi1 = 0.93125 / 0.7625 and phi2 = (0.25 + 0.3 phi1) / 0.8.
  const auto phi = ExactFlux();
  const double phi1 = 0.93125 / 0.7625;
  EXPECT_NEAR(phi[0], 2.0, 1.0e-14);
  EXPECT_NEAR(phi[1], phi1, 1.0e-14);
  EXPECT_NEAR(phi[2], (0.25 + 0.3 * phi1) / 0.8, 1.0e-14);
}

TEST_P(GroupsetSolveFluxViewsTest, ReproducesTheInfiniteMediumFlux)
{
  const auto& [ranges, method] = GetParam();
  auto problem = MakeProblem(ranges, method, true);
  for (std::size_t gs = 0; gs < ranges.size(); ++gs)
    ASSERT_GT(GlobalDelayedUnknowns(*problem, gs), 0U) << "groupset " << gs;

  // The last block's delayed unknowns fit in the tail, on every rank.
  const auto& last_groupset = problem->GetGroupsets().back();
  ASSERT_EQ(problem->GetFluxBlock(last_groupset), ranges.size() - 1);
  const auto local_delayed = last_groupset.angle_agg->GetNumDelayedAngularDOFs().first;
  EXPECT_EQ(problem->GetPhiOldLocal().GetTailCapacity(), local_delayed);
  EXPECT_EQ(problem->GetPhiNewLocal().GetTailCapacity(), local_delayed);

  const auto& q_moments = problem->GetQMomentsLocal();
  const std::vector<double> q_before(q_moments.Data().begin(), q_moments.Data().end());

  SolveGroupsets(*problem, false);
  EXPECT_LT(MaxRelativeFluxError(*problem), kAcceptance);
  for (std::size_t gs = 0; gs < ranges.size(); ++gs)
    EXPECT_EQ(problem->GetWGSContext(static_cast<int>(gs)).last_solve.status,
              IterationStatus::CONVERGED)
      << "groupset " << gs;

  // Solving again, with the solution as the initial guess, keeps it.
  SolveGroupsets(*problem, false);
  EXPECT_LT(MaxRelativeFluxError(*problem), kAcceptance);

  // So does solving again after the workspace is freed and rebuilt, from a zero flux.
  for (std::size_t gs = 0; gs < ranges.size(); ++gs)
    problem->GetWGSSolver(gs)->ReleaseWorkspace();
  problem->GetPhiOldLocal().Fill(0.0);
  problem->GetPhiNewLocal().Fill(0.0);
  SolveGroupsets(*problem, true);
  EXPECT_LT(MaxRelativeFluxError(*problem), kAcceptance);

  const auto q_after = q_moments.Data();
  ASSERT_EQ(q_after.size(), q_before.size());
  EXPECT_TRUE(std::equal(q_after.begin(), q_after.end(), q_before.begin()));
}

INSTANTIATE_TEST_SUITE_P(
  InnerMethods,
  GroupsetSolveFluxViewsTest,
  testing::Combine(testing::Values(GroupsetRanges{{0, 2}}, GroupsetRanges{{0, 0}, {1, 2}}),
                   testing::Values("petsc_gmres", "petsc_bicgstab", "petsc_richardson")),
  [](const testing::TestParamInfo<SolveCase>& info)
  {
    return std::get<1>(info.param) + "_" + std::to_string(std::get<0>(info.param).size()) +
           "_groupsets";
  });

/*
 * A solve starts from the current flux unless that flux is zero. The solution vector views the
 * flux storage, which the problem writes outside PETSc, so the solver must not use a norm that
 * PETSc cached before such a write. Solve with a zero source, which leaves a zero flux, then write
 * the converged solution of the nonzero-source problem, restore the source, and solve again: the
 * solve must start from that solution and take no iterations.
 */
namespace
{

class GroupsetSolveInitialGuessTest : public testing::TestWithParam<bool>
{
};

} // namespace

TEST_P(GroupsetSolveInitialGuessTest, StartsFromAFluxWrittenOutsidePETSc)
{
  const bool reflecting = GetParam();
  auto problem = MakeProblem({{0, 2}}, "petsc_gmres", reflecting);
  auto& angle_agg = *problem->GetGroupsets().front().angle_agg;
  EXPECT_EQ(GlobalDelayedUnknowns(*problem, 0) > 0, reflecting);
  auto solver = problem->GetWGSSolver(0);
  const auto& context = problem->GetWGSContext(0);

  solver->Setup();
  solver->Solve();
  ASSERT_EQ(context.last_solve.status, IterationStatus::CONVERGED);
  ASSERT_GT(context.last_solve.num_iterations, 1U);
  if (reflecting)
  {
    ASSERT_LT(MaxRelativeFluxError(*problem), kAcceptance);
  }
  const auto phi = std::as_const(problem->GetPhiOldLocal()).Data();
  const std::vector<double> solution(phi.begin(), phi.end());
  const auto delayed_solution = angle_agg.GetOldDelayedAngularDOFsAsSTLVector();

  problem->SetVolumetricSources({MakeSource(0.0)}, true);
  problem->GetPhiOldLocal().Fill(0.0);
  angle_agg.SetOldDelayedAngularDOFsFromSTLVector(
    std::vector<double>(delayed_solution.size(), 0.0));
  solver->Solve();
  ASSERT_EQ(context.last_solve.status, IterationStatus::CONVERGED);
  const auto zero_phi = std::as_const(problem->GetPhiOldLocal()).Data();
  ASSERT_TRUE(std::all_of(zero_phi.begin(), zero_phi.end(), [](double v) { return v == 0.0; }));

  problem->SetVolumetricSources({MakeSource(1.0)}, true);
  problem->GetPhiOldLocal().Assign(solution);
  angle_agg.SetOldDelayedAngularDOFsFromSTLVector(delayed_solution);
  solver->Solve();
  EXPECT_EQ(context.last_solve.status, IterationStatus::CONVERGED);
  EXPECT_EQ(context.last_solve.num_iterations, 0U);
  if (reflecting)
  {
    EXPECT_LT(MaxRelativeFluxError(*problem), kAcceptance);
  }
}

INSTANTIATE_TEST_SUITE_P(Boundaries,
                         GroupsetSolveInitialGuessTest,
                         testing::Bool(),
                         [](const testing::TestParamInfo<bool>& info)
                         {
                           return info.param ? std::string("reflecting") : std::string("vacuum");
                         });
