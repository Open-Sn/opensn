// SPDX-FileCopyrightText: 2026 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/solvers/solver_scheme.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/discrete_ordinates_problem.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/iterative_methods/ags_linear_solver.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/iterative_methods/classic_richardson.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/iterative_methods/sweep_wgs_context.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/iterative_methods/wgs_linear_solver.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep/angle_set/angle_set.h"
#include "modules/linear_boltzmann_solvers/lbs_problem/groupset/lbs_groupset.h"
#include "framework/utils/error.h"
#include <algorithm>

namespace opensn
{

namespace
{

void
UpdateSweepWorkLimits(SolverScheme& scheme, LBSGroupset& groupset)
{
  scheme.max_groupset_size = std::max(scheme.max_groupset_size, groupset.GetNumGroups());

  for (auto& angleset : *groupset.angle_agg)
  {
    const auto& spds = angleset->GetSPDS();
    const auto& levelized_spls = spds.GetLevelizedLocalSubgrid();
    for (const auto& level : levelized_spls)
      scheme.max_level_size = std::max(scheme.max_level_size, level.size());

    scheme.max_angleset_size =
      std::max(scheme.max_angleset_size, angleset->GetAngleIndices().size());
  }
}

std::shared_ptr<WGSContext>
MakeWGSContext(DiscreteOrdinatesProblem& problem,
               LBSGroupset& groupset,
               const SetSourceFunction& set_source_function,
               const SweepChunkFactory& sweep_chunk_factory)
{
  const auto& options = problem.GetOptions();
  return std::make_shared<SweepWGSContext>(problem,
                                           groupset,
                                           set_source_function,
                                           APPLY_WGS_SCATTER_SOURCES | APPLY_WGS_FISSION_SOURCES,
                                           APPLY_FIXED_SOURCES | APPLY_AGS_SCATTER_SOURCES |
                                             APPLY_AGS_FISSION_SOURCES,
                                           options.verbose_inner_iterations,
                                           sweep_chunk_factory(groupset));
}

std::shared_ptr<LinearSolver>
MakeWGSSolver(LBSGroupset& groupset,
              const std::shared_ptr<WGSContext>& wgs_context,
              bool verbose_inner_iterations)
{
  if (groupset.iterative_method == LinearSystemSolver::IterativeMethod::CLASSIC_RICHARDSON)
    return std::make_shared<ClassicRichardson>(wgs_context, verbose_inner_iterations);

  return std::make_shared<WGSLinearSolver>(wgs_context);
}

std::shared_ptr<AGSLinearSolver>
MakeAGSSolver(DiscreteOrdinatesProblem& problem,
              const std::vector<std::shared_ptr<LinearSolver>>& wgs_solvers)
{
  const auto& options = problem.GetOptions();
  auto ags_solver = std::make_shared<AGSLinearSolver>(problem, wgs_solvers);
  if (problem.GetNumGroupsets() == 1)
  {
    ags_solver->SetMaxIterations(1);
    ags_solver->SetVerbosity(false);
  }
  else
  {
    ags_solver->SetMaxIterations(options.max_ags_iterations);
    ags_solver->SetVerbosity(options.verbose_inner_iterations);
  }
  ags_solver->SetTolerance(options.ags_tolerance);

  return ags_solver;
}

} // namespace

SolverScheme
BuildSolvers(DiscreteOrdinatesProblem& problem,
             const SetSourceFunction& set_source_function,
             const SweepChunkFactory& sweep_chunk_factory)
{
  SolverScheme result;
  const auto& options = problem.GetOptions();

  for (size_t gsid = 0; gsid < problem.GetNumGroupsets(); ++gsid)
  {
    auto& groupset = problem.GetGroupset(gsid);
    UpdateSweepWorkLimits(result, groupset);
    result.wgs_contexts.push_back(
      MakeWGSContext(problem, groupset, set_source_function, sweep_chunk_factory));
  }

  OpenSnLogicalErrorIf(result.wgs_contexts.empty(),
                       problem.GetName() + ": Cannot initialize WGS solvers before WGS contexts.");
  OpenSnLogicalErrorIf(result.wgs_contexts.size() != problem.GetNumGroupsets(),
                       problem.GetName() + ": WGS context count does not match groupset count.");

  for (size_t gsid = 0; gsid < problem.GetNumGroupsets(); ++gsid)
  {
    auto& groupset = problem.GetGroupset(gsid);
    const auto& wgs_context = result.wgs_contexts[gsid];
    OpenSnLogicalErrorIf(not wgs_context, problem.GetName() + ": Null WGS context.");
    result.wgs_solvers.push_back(
      MakeWGSSolver(groupset, wgs_context, options.verbose_inner_iterations));
  }

  result.ags_solver = MakeAGSSolver(problem, result.wgs_solvers);
  return result;
}

} // namespace opensn
