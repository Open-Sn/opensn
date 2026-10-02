// SPDX-FileCopyrightText: 2026 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "test/unit/common/mesh_builders.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/acceleration/cmfd_acceleration.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/discrete_ordinates_problem.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/solvers/pi_keigen_solver.h"
#include "framework/materials/multi_group_xs/multi_group_xs.h"
#include "framework/math/quadratures/angular/product_quadrature.h"
#include "framework/parameters/parameter_block.h"
#include "framework/runtime.h"
#include <gtest/gtest.h>
#include <chrono>
#include <cmath>
#include <cstddef>
#include <filesystem>
#include <fstream>
#include <memory>
#include <string>

using namespace opensn;

namespace
{

/// Loads a one-group fissile material (sigma_t = 1, c = 0.9, nu sigma_f = 0.12).
std::shared_ptr<MultiGroupXS>
MakeFissileXS()
{
  const auto now = std::chrono::steady_clock::now().time_since_epoch().count();
  const auto path = std::filesystem::temp_directory_path() /
                    ("opensn_cmfd_repeated_execute_r" + std::to_string(mpi_comm.rank()) + "_" +
                     std::to_string(now) + ".xs");
  {
    std::ofstream out(path);
    out << "NUM_GROUPS 1\nNUM_MOMENTS 1\n"
        << "SIGMA_T_BEGIN\n0 1.0\nSIGMA_T_END\n"
        << "SIGMA_F_BEGIN\n0 0.05\nSIGMA_F_END\n"
        << "NU_PROMPT_BEGIN\n0 2.4\nNU_PROMPT_END\n"
        << "CHI_PROMPT_BEGIN\n0 1.0\nCHI_PROMPT_END\n"
        << "TRANSFER_MOMENTS_BEGIN\nM_GFROM_GTO_VAL 0 0 0 0.9\nTRANSFER_MOMENTS_END\n";
  }
  auto xs = std::make_shared<MultiGroupXS>(MultiGroupXS::LoadFromOpenSn(path.string()));
  std::filesystem::remove(path);
  return xs;
}

/// A 2D XY problem with one classic-Richardson groupset of `l_max_its` iterations.
std::shared_ptr<DiscreteOrdinatesProblem>
MakeProblem(const int l_max_its,
            const bool adjoint = false,
            const unsigned int scattering_order = 0,
            const std::string& write_restart_stem = "",
            const std::string& read_restart_stem = "")
{
  auto grid = BuildSquareMesh(10.0, 4, 0.0);
  grid->SetUniformBlockID(0);

  ParameterBlock groupset("0");
  groupset.AddParameter(ParameterBlock("groups_from_to", std::vector<int>{0, 0}));
  groupset.AddParameter("angular_quadrature",
                        std::shared_ptr<AngularQuadrature>(
                          std::make_shared<GLCProductQuadrature2DXY>(2, 4, scattering_order)));
  groupset.AddParameter("inner_linear_method", std::string("classic_richardson"));
  groupset.AddParameter("l_max_its", l_max_its);
  groupset.AddParameter("l_abs_tol", 1.0e-14);
  ParameterBlock groupsets("groupsets");
  groupsets.ChangeToArray();
  groupsets.AddParameter(groupset);

  ParameterBlock xs_entry("0");
  xs_entry.AddParameter(ParameterBlock("block_ids", std::vector<int>{0}));
  xs_entry.AddParameter("xs", MakeFissileXS());
  ParameterBlock xs_map("xs_map");
  xs_map.ChangeToArray();
  xs_map.AddParameter(xs_entry);

  ParameterBlock params;
  params.AddParameter("mesh", grid);
  params.AddParameter("num_groups", 1);
  params.AddParameter(groupsets);
  params.AddParameter(xs_map);
  ParameterBlock options("options");
  options.AddParameter("adjoint", adjoint);
  options.AddParameter("save_angular_flux", true);
  if (not write_restart_stem.empty())
    options.AddParameter("write_restart_path", write_restart_stem);
  if (not read_restart_stem.empty())
    options.AddParameter("read_restart_path", read_restart_stem);
  params.AddParameter(options);
  return DiscreteOrdinatesProblem::Create(params);
}

/// Executes the solver and returns the sweeps it did per power iteration.
double
SweepsPerOuter(PowerIterationKEigenSolver& solver)
{
  const auto sweeps_before = solver.GetNumSweeps();
  solver.Execute();
  return static_cast<double>(solver.GetNumSweeps() - sweeps_before) /
         solver.GetNumPowerIterations();
}

} // namespace

// CMFD is inactive for the whole solve (inactive_iterations > max_iters), so every transport solve
// uses the caller's within-group settings: classic Richardson with an unreachable tolerance does
// exactly l_max_its sweeps per power iteration. Changing the settings between two solves of the
// same solver must take effect in the second solve, and both solves must leave them unchanged.
TEST(CMFDRepeatedExecuteTest, UsesCurrentWGSSettings)
{
  auto problem = MakeProblem(2);

  ParameterBlock cmfd_params;
  cmfd_params.AddParameter("problem", std::shared_ptr<Problem>(problem));
  cmfd_params.AddParameter("coarse_mesh", std::string("identity"));
  cmfd_params.AddParameter("inactive_iterations", 1000);
  auto cmfd = CMFDAcceleration::Create(cmfd_params);

  ParameterBlock pi_params;
  pi_params.AddParameter("problem", std::shared_ptr<Problem>(problem));
  pi_params.AddParameter("acceleration",
                         std::shared_ptr<DiscreteOrdinatesKEigenAcceleration>(cmfd));
  pi_params.AddParameter("max_iters", 3);
  pi_params.AddParameter("k_tol", 1.0e-14);
  auto solver = PowerIterationKEigenSolver::Create(pi_params);
  solver->Initialize();

  auto& groupset = problem->GetGroupset(0);
  EXPECT_DOUBLE_EQ(SweepsPerOuter(*solver), 2.0);
  EXPECT_EQ(groupset.max_iterations, 2U);

  problem->SetWGSSolveSettings(static_cast<std::size_t>(groupset.id), 3, 1.0e-14);
  EXPECT_DOUBLE_EQ(SweepsPerOuter(*solver), 3.0);
  EXPECT_EQ(groupset.max_iterations, 3U);
}

// The first eigenvalue change of a solve is measured from the solver's current eigenvalue. After a
// converged solve, a second Execute() starts at the converged k: its first two changes are below
// k_tol, so it converges in exactly two power iterations (the minimum for two consecutive small
// changes). Measuring the first change from 1.0 instead would take three, and could count an
// unrelated small first change toward convergence.
TEST(PIKEigenRepeatedExecuteTest, FirstChangeIsMeasuredFromCurrentEigenvalue)
{
  auto problem = MakeProblem(2);

  ParameterBlock pi_params;
  pi_params.AddParameter("problem", std::shared_ptr<Problem>(problem));
  pi_params.AddParameter("max_iters", 1000);
  pi_params.AddParameter("k_tol", 1.0e-8);
  auto solver = PowerIterationKEigenSolver::Create(pi_params);
  solver->Initialize();

  solver->Execute();
  const double k_first = solver->GetEigenvalue();
  ASSERT_LT(solver->GetNumPowerIterations(), 1000U);
  // Otherwise a change measured from 1.0 would also be small and the test would not discriminate.
  ASSERT_GT(std::fabs(k_first - 1.0), 1.0e-3);

  solver->Execute();
  EXPECT_EQ(solver->GetNumPowerIterations(), 2U);
  EXPECT_NEAR(solver->GetEigenvalue(), k_first, 1.0e-8 * k_first);
}

TEST(PIKEigenRepeatedExecuteTest, PhysicalAdjointRestartIsNormalizedForSolver)
{
  auto now = std::chrono::steady_clock::now().time_since_epoch().count();
  mpi_comm.broadcast(now, 0);
  const auto restart_stem = (std::filesystem::temp_directory_path() /
                             ("opensn_adjoint_physical_restart_" + std::to_string(now) + "_"))
                              .string();
  const auto restart_path = restart_stem + std::to_string(mpi_comm.rank()) + ".restart.h5";

  auto problem = MakeProblem(20, true, 1, restart_stem);
  ParameterBlock pi_params;
  pi_params.AddParameter("problem", std::shared_ptr<Problem>(problem));
  pi_params.AddParameter("max_iters", 1000);
  pi_params.AddParameter("k_tol", 1.0e-10);
  auto solver = PowerIterationKEigenSolver::Create(pi_params);
  solver->Initialize();
  solver->Execute();
  const auto physical_phi = problem->GetPhiNewLocal();
  ASSERT_TRUE(problem->IsAdjointSolutionReoriented());
  ASSERT_TRUE(problem->WriteRestartData());

  auto restarted_problem = MakeProblem(20, true, 1, "", restart_stem);
  ParameterBlock restarted_pi_params;
  restarted_pi_params.AddParameter("problem", std::shared_ptr<Problem>(restarted_problem));
  restarted_pi_params.AddParameter("max_iters", 1000);
  restarted_pi_params.AddParameter("k_tol", 1.0e-10);
  auto restarted_solver = PowerIterationKEigenSolver::Create(restarted_pi_params);
  restarted_solver->Initialize();

  EXPECT_FALSE(restarted_problem->IsAdjointSolutionReoriented());
  restarted_problem->ReorientAdjointSolution();
  ASSERT_EQ(restarted_problem->GetPhiNewLocal().size(), physical_phi.size());
  for (size_t i = 0; i < physical_phi.size(); ++i)
    EXPECT_NEAR(restarted_problem->GetPhiNewLocal()[i], physical_phi[i], 1.0e-12);

  restarted_problem->RestoreAdjointSolutionOrientation();
  EXPECT_FALSE(restarted_problem->IsAdjointSolutionReoriented());

  // Independently exercise Y_{ell,m}(-Omega) = (-1)^ell Y_{ell,m}(Omega) with known
  // moment values. This prevents the restart round-trip above from passing if the
  // orientation transform is accidentally an identity operation.
  const auto& moment_map =
    restarted_problem->GetGroupset(0).quadrature->GetMomentToHarmonicsIndexMap();
  ASSERT_GT(moment_map.size(), 1U);
  ASSERT_EQ(moment_map[0].ell, 0U);
  ASSERT_EQ(moment_map[1].ell, 1U);
  const auto& local_cells = restarted_problem->GetGrid()->GetLocalCells();
  if (not local_cells.empty())
  {
    const auto& cell = local_cells.front();
    const auto& sdm = restarted_problem->GetSpatialDiscretization();
    const auto& uk_man = restarted_problem->GetUnknownManager();
    const auto even_dof = sdm.MapDOFLocal(cell, 0, uk_man, 0, 0);
    const auto odd_dof = sdm.MapDOFLocal(cell, 0, uk_man, 1, 0);
    auto& phi_new = restarted_problem->GetPhiNewLocal();
    auto& phi_old = restarted_problem->GetPhiOldLocal();
    phi_new[even_dof] = 2.5;
    phi_new[odd_dof] = 3.25;
    phi_old[even_dof] = -1.5;
    phi_old[odd_dof] = -4.75;

    restarted_problem->ReorientAdjointSolution();
    EXPECT_DOUBLE_EQ(phi_new[even_dof], 2.5);
    EXPECT_DOUBLE_EQ(phi_new[odd_dof], -3.25);
    EXPECT_DOUBLE_EQ(phi_old[even_dof], -1.5);
    EXPECT_DOUBLE_EQ(phi_old[odd_dof], 4.75);
  }

  std::filesystem::remove(restart_path);
}
