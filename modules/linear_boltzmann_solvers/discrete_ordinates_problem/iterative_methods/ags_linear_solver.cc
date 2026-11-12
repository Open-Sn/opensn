// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/iterative_methods/ags_linear_solver.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/iterative_methods/convergence.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/iterative_methods/sweep_wgs_context.h"
#include "modules/linear_boltzmann_solvers/lbs_problem/iterative_methods/iteration_logging.h"
#include "modules/linear_boltzmann_solvers/lbs_problem/lbs_problem.h"
#include "framework/mesh/mesh/mesh.h"
#include "framework/logging/log.h"
#include "framework/utils/caliper_scopes.h"
#include "framework/runtime.h"
#include "caliper/cali.h"
#include <algorithm>
#include <cmath>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

namespace opensn
{

namespace
{

/// Throws on all ranks if the scalar flux contains a non-finite value after a groupset solve.
void
ValidateFiniteScalarFlux(const LBSProblem& lbs_problem, const std::shared_ptr<LinearSolver>& solver)
{
  const auto& phi = lbs_problem.GetPhiNewLocal();
  const int local_nonfinite =
    std::all_of(phi.begin(), phi.end(), [](const double value) { return std::isfinite(value); })
      ? 0
      : 1;
  int global_nonfinite = 0;
  mpi_comm.all_reduce(local_nonfinite, global_nonfinite, mpi::op::max<int>());
  if (global_nonfinite == 0)
    return;

  std::string groups;
  if (const auto context = std::dynamic_pointer_cast<WGSContext>(solver->GetContext()))
    groups = " for groups [" + std::to_string(context->groupset.first_group) + "-" +
             std::to_string(context->groupset.last_group) + "]";
  throw std::runtime_error(
    "The within-group solve" + groups +
    " produced a non-finite scalar flux. Check the cross sections, sources, and boundary "
    "conditions for invalid values or regions without interaction.");
}

/**
 * Returns whether one pass over the groupsets, in the order of wgs_solvers, solves the coupling
 * between them exactly. A groupset receives sources from the others through scattering and, when
 * its right-hand side applies across-groupset fission, through fission. The pass is exact when no
 * groupset receives a source from one solved after it: no material scatters from a later groupset
 * into an earlier one, in the problem's mode (transposed in adjoint mode), and either no material
 * is fissionable or no groupset applies across-groupset fission. Continuous slowing down couples
 * groups through the sweep, so a CSDA problem always iterates. Every rank holds the cross
 * sections of every block, so every rank gets the same answer.
 */
bool
OnePassIsExact(const LBSProblem& lbs_problem,
               const std::vector<std::shared_ptr<LinearSolver>>& wgs_solvers)
{
  const auto& options = lbs_problem.GetOptions();
  if (options.csda_enabled)
    return false;

  std::vector<size_t> groupset_of(lbs_problem.GetNumGroups(), wgs_solvers.size());
  bool across_groupset_fission = false;
  for (size_t k = 0; k < wgs_solvers.size(); ++k)
  {
    const auto context = std::dynamic_pointer_cast<WGSContext>(wgs_solvers[k]->GetContext());
    if (not context)
      return false;
    for (auto g = context->groupset.first_group; g <= context->groupset.last_group; ++g)
      groupset_of[g] = k;
    across_groupset_fission =
      across_groupset_fission or (context->GetRHSSourceScope() & APPLY_AGS_FISSION_SOURCES);
  }
  if (std::find(groupset_of.begin(), groupset_of.end(), wgs_solvers.size()) != groupset_of.end())
    return false;

  for (const auto& [block_id, xs] : lbs_problem.GetBlockID2XSMap())
  {
    if (across_groupset_fission and xs->IsFissionable())
      return false;
    for (const auto& transfer_matrix : xs->GetTransferMatrices(options.adjoint))
      for (size_t g = 0; g < transfer_matrix.GetNumRows(); ++g)
        for (const auto& [row, source_group, value] : transfer_matrix.Row(g))
          if (value != 0.0 and groupset_of[source_group] > groupset_of[g])
            return false;
  }
  return true;
}

/// L2 norm of a distributed vector over all ranks.
double
ComputeGlobalL2Norm(const std::vector<double>& values)
{
  double local_sum = 0.0;
  for (const double value : values)
    local_sum += value * value;
  double global_sum = 0.0;
  mpi_comm.all_reduce<double>(local_sum, global_sum, mpi::op::sum<double>());
  return std::sqrt(global_sum);
}

} // namespace

void
AGSLinearSolver::Solve()
{
  CaliperPhaseScope cali_solve_phase("Solve", CaliperSolvePhaseDepth());
  CaliperRegionScope cali_ags("AGS", CaliperAGSScopeDepth());
  const bool is_multi_groupset = wgs_solvers_.size() > 1;
  const bool print_ags_stats = verbose_ and is_multi_groupset;
  last_solve_ = {};

  phi_old_ = lbs_problem_.GetPhiOldLocal();

  // Save qmoms to be restored after each iteration. This is necessary for multiple ags iterations
  // to function and for keigen-value problems
  const auto saved_qmoms = lbs_problem_.GetQMomentsLocal();

  // Without coupling to earlier groupsets, a second pass would only continue the groupset solves
  // that already converged, so the iteration stops after the first pass in which they all do.
  const bool one_pass_exact = is_multi_groupset and OnePassIsExact(lbs_problem_, wgs_solvers_);

  ConvergenceRateEstimator rate_estimator;
  double last_error = 0.0;
  bool converged = false;
  bool failed = false;
  unsigned int num_iterations = 0;
  const bool reports_as_iteration = is_multi_groupset;
  CALI_CXX_MARK_LOOP_BEGIN(ags_iteration, "AGSIteration");
  for (unsigned int iter = 0; iter < max_iterations_; ++iter)
  {
    CALI_CXX_MARK_LOOP_ITERATION(ags_iteration, iter);

    num_iterations = iter + 1;
    for (auto& solver : wgs_solvers_)
    {
      solver->Setup();
      solver->Solve();
      ValidateFiniteScalarFlux(lbs_problem_, solver);
    }

    std::vector<IterationSummary> wgs_summaries;
    wgs_summaries.reserve(wgs_solvers_.size());
    for (size_t gsid = 0; gsid < wgs_solvers_.size(); ++gsid)
    {
      auto wgs_context = std::dynamic_pointer_cast<WGSContext>(wgs_solvers_[gsid]->GetContext());
      if (wgs_context and HasIterationStatus(wgs_context->last_solve))
        wgs_summaries.emplace_back(wgs_context->last_solve);
    }
    const auto wgs_status = MostSevereIterationStatus(wgs_summaries);

    std::stringstream iter_stats;
    iter_stats << program_timer.GetTimeString() << " AGS iteration = " << iter;

    if (lbs_problem_.GetOptions().ags_pointwise_convergence)
    {
      double pw_change = ComputePointwisePhiChange(lbs_problem_, phi_old_);
      last_error = pw_change;
      const double rho = rate_estimator.Update(pw_change);

      AppendNumericField(iter_stats, "pw_change", pw_change, Scientific(6));
      AppendNumericField(iter_stats, "rho_est", rho, Fixed(4));

      if (rate_estimator.Converged(tolerance_))
        converged = true;
    }
    else
    {
      double norm = ComputeL2PhiChange(lbs_problem_, phi_old_);
      last_error = norm;
      const double rho = rate_estimator.Update(norm);

      AppendNumericField(iter_stats, "l2_change", norm, Scientific(6));
      AppendNumericField(iter_stats, "rho_est", rho, Fixed(4));

      // The L2 change is absolute, so its round-off level scales with the flux norm.
      const double roundoff_change = ConvergenceRateEstimator::kRoundoffRelativeChange *
                                     ComputeGlobalL2Norm(lbs_problem_.GetPhiNewLocal());
      if (rate_estimator.Converged(tolerance_, roundoff_change))
        converged = true;
    }

    if (one_pass_exact and wgs_summaries.size() == wgs_solvers_.size() and
        wgs_status == IterationStatus::CONVERGED)
      converged = true;

    const auto iteration_status = IterationStatusFromSolve(converged, false, wgs_status);
    if (iteration_status != IterationStatus::NONE)
      iter_stats << ", status = " << IterationStatusName(iteration_status);
    iter_stats << FormatNestedStatusCounts("WGS", wgs_summaries);

    if (print_ags_stats)
      log.Log() << no_wrap << iter_stats.str();

    // Restore qmoms
    lbs_problem_.SetQMomentsFrom(saved_qmoms);

    // A within-group solve that stops at its iteration limit is an inexact inner solve; the
    // across-groupset iteration continues and its own convergence test decides.
    if (wgs_status == IterationStatus::FAILED)
    {
      failed = true;
      break;
    }
    else if (converged)
      break;
    else
      phi_old_ = lbs_problem_.GetPhiNewLocal();
  }
  CALI_CXX_MARK_LOOP_END(ags_iteration);

  const auto inner_status = failed ? IterationStatus::FAILED : IterationStatus::NONE;
  last_solve_.num_iterations = reports_as_iteration ? num_iterations : 0;
  last_solve_.status = reports_as_iteration
                         ? IterationStatusFromSolve(converged, true, inner_status)
                         : IterationStatus::NONE;
  last_solve_.metric_name =
    lbs_problem_.GetOptions().ags_pointwise_convergence ? "pw_change" : "l2_change";
  last_solve_.metric_value = last_error;
  if (print_ags_stats and (last_solve_.status == IterationStatus::FAILED or
                           last_solve_.status == IterationStatus::LIMIT))
    log.Log() << no_wrap << program_timer.GetTimeString() << " "
              << FormatIterationSummary("AGS", last_solve_);

  for (const auto& solver : wgs_solvers_)
  {
    auto sweep_context = std::dynamic_pointer_cast<SweepWGSContext>(solver->GetContext());
    if (not sweep_context or not sweep_context->log_info)
      continue;

    const auto stats = sweep_context->GetAccumulatedSweepStats();
    if (stats.num_sweeps == 0)
      continue;

    const double avg_sweep_time = stats.total_sweep_time / static_cast<double>(stats.num_sweeps);
    const auto& groupset = sweep_context->groupset;
    const size_t num_angles = groupset.quadrature->GetNumAngles();
    const size_t num_unknowns =
      lbs_problem_.GetGlobalNodeCount() * num_angles * groupset.GetNumGroups();
    const auto num_delayed_psi_info = groupset.angle_agg->GetNumDelayedAngularDOFs();
    const size_t num_delayed_unknowns = num_delayed_psi_info.second;
    const double delayed_unknown_percent =
      (num_unknowns > 0)
        ? static_cast<double>(num_delayed_unknowns) * 100.0 / static_cast<double>(num_unknowns)
        : 0.0;
    double max_avg_sweep_time = 0.0;
    opensn::mpi_comm.all_reduce(&avg_sweep_time, 1, &max_avg_sweep_time, mpi::op::max<double>());

    const double sweep_time_per_unknown =
      num_unknowns > 0 ? max_avg_sweep_time * 1.0e9 / static_cast<double>(num_unknowns) : 0.0;
    const std::string label = "WGS groups [" + std::to_string(groupset.first_group) + "-" +
                              std::to_string(groupset.last_group) + "]";
    std::stringstream sweep_timing;
    sweep_timing << label;
    AppendNumericField(sweep_timing, "avg_sweep_time", max_avg_sweep_time, Scientific(6), false);
    sweep_timing << " s";
    AppendNumericField(
      sweep_timing, "sweep_time_per_unknown", sweep_time_per_unknown, Scientific(6));
    sweep_timing << " ns";
    log.Log() << no_wrap << program_timer.GetTimeString() << " " << sweep_timing.str();

    std::stringstream sweep_work;
    sweep_work << label;
    AppendNumericField(sweep_work, "unknowns", num_unknowns, false);
    AppendNumericField(sweep_work, "lagged_unknowns", num_delayed_unknowns);
    AppendNumericField(sweep_work, "lagged_pct", delayed_unknown_percent, Fixed(2));
    log.Log() << no_wrap << program_timer.GetTimeString() << " " << sweep_work.str();

    sweep_context->ResetAccumulatedSweepStats();
  }
}

} // namespace opensn
