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
#include <numeric>
#include <optional>
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
 * Checks whether the across-groupset iteration must iterate. Returns why one pass over the
 * groupsets, in the solve order of wgs_solvers, does not solve the coupling between them exactly,
 * or nothing when it does. A groupset receives sources from the
 * others through scattering and, when its right-hand side applies across-groupset fission, through
 * fission. The pass is exact when no groupset receives a source from one solved after it: no
 * material scatters from a later groupset into an earlier one, in the problem's mode (transposed in
 * adjoint mode), and either no material is fissionable or no groupset applies across-groupset
 * fission. A CSDA problem always iterates. Its reason suggests no other groupsets, since each
 * contiguous block of charged-particle groups must stay in one groupset.
 *
 * When scattering blocks one pass, the reason names the first such transfer found and how to remove
 * all of them: the reverse order, when that is exact with the same groupsets, or otherwise the
 * group ranges that must each share a groupset in this order. A transfer from group g' into group g
 * against the solve order, upward in energy in ascending order and downward in descending order,
 * is solved in one pass only when g, g', and the groups between them share a groupset. Every rank
 * holds the cross sections of every block in the same order, so every rank gets the same answer.
 */
std::optional<std::string>
CheckAGSIterates(const LBSProblem& lbs_problem,
                 const std::vector<std::shared_ptr<LinearSolver>>& wgs_solvers)
{
  const auto& options = lbs_problem.GetOptions();
  if (options.csda_enabled)
    return "CSDA is enabled.";

  const size_t num_solvers = wgs_solvers.size();
  std::vector<size_t> solve_index(lbs_problem.GetNumGroups(), num_solvers);
  bool across_groupset_fission = false;
  for (size_t k = 0; k < num_solvers; ++k)
  {
    const auto context = std::dynamic_pointer_cast<WGSContext>(wgs_solvers[k]->GetContext());
    if (not context)
      return "a groupset has no within-group solver.";
    for (auto g = context->groupset.first_group; g <= context->groupset.last_group; ++g)
      solve_index[g] = k;
    across_groupset_fission =
      across_groupset_fission or (context->GetRHSSourceScope() & APPLY_AGS_FISSION_SOURCES);
  }
  if (std::find(solve_index.begin(), solve_index.end(), num_solvers) != solve_index.end())
    return "the groupsets do not cover every group.";

  const auto& xs_map = lbs_problem.GetBlockID2XSMap();
  if (across_groupset_fission)
    for (const auto& [_, xs] : xs_map)
      if (xs->IsFissionable())
        return "fission couples the groupsets.";

  // Cross sections may have more groups than the problem simulates. Only transfers between
  // simulated groups couple the groupsets.
  const size_t num_groups = solve_index.size();
  const bool ascending = solve_index.front() < solve_index.back();
  std::optional<std::string> first_blocking;
  bool reverse_is_exact = true;
  // The highest group that must share a groupset with each group, from the transfers against the
  // solve order. A transfer against the order joins two different groups, so reach[g] = g means
  // that group g has none.
  std::vector<size_t> reach(num_groups);
  std::iota(reach.begin(), reach.end(), size_t{0});
  for (const auto& [_, xs] : xs_map)
    for (const auto& transfer_matrix : xs->GetTransferMatrices(options.adjoint))
      for (size_t g = 0; g < std::min<size_t>(transfer_matrix.GetNumRows(), num_groups); ++g)
        for (const auto& [row, source_group, value] : transfer_matrix.Row(g))
        {
          if (value == 0.0 or source_group >= num_groups)
            continue;
          if (solve_index[source_group] > solve_index[g] and not first_blocking)
            first_blocking = "group " + std::to_string(source_group) + " scatters into group " +
                             std::to_string(g) + " in an earlier groupset" +
                             (options.adjoint ? " (adjoint)." : ".");
          if (solve_index[source_group] < solve_index[g])
            reverse_is_exact = false;
          if (ascending ? source_group > g : source_group < g)
          {
            const size_t low = std::min<size_t>(source_group, g);
            reach[low] = std::max<size_t>({reach[low], source_group, g});
          }
        }
  if (not first_blocking)
    return std::nullopt;

  if (reverse_is_exact)
  {
    // The reverse of the current order, which "auto" gives when it is the mode's default.
    const bool reverse_is_descending = ascending;
    const bool auto_is_descending = options.adjoint;
    const std::string order = reverse_is_descending == auto_is_descending
                                ? "auto"
                                : (reverse_is_descending ? "descending" : "ascending");
    return *first_blocking + " Set ags_groupset_order = \"" + order + "\" for a single pass.";
  }

  // Merge the ranges that must share a groupset.
  std::vector<std::pair<size_t, size_t>> shared_ranges;
  for (size_t g = 0; g < num_groups; ++g)
  {
    if (reach[g] == g)
      continue;
    if (not shared_ranges.empty() and g <= shared_ranges.back().second)
      shared_ranges.back().second = std::max(shared_ranges.back().second, reach[g]);
    else
      shared_ranges.emplace_back(g, reach[g]);
  }
  if (shared_ranges.size() == 1 and shared_ranges.front().first == 0 and
      shared_ranges.front().second + 1 == num_groups)
    return *first_blocking + " Use one groupset for a single pass.";
  std::string ranges;
  for (size_t i = 0; i < shared_ranges.size(); ++i)
  {
    if (i > 0)
      ranges += i + 1 == shared_ranges.size() ? " and " : ", ";
    ranges +=
      std::to_string(shared_ranges[i].first) + "-" + std::to_string(shared_ranges[i].second);
  }
  return *first_blocking + " Put groups " + ranges +
         (shared_ranges.size() == 1 ? " in one groupset" : " each in one groupset") +
         " for a single pass.";
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
  const auto& options = lbs_problem_.GetOptions();
  const bool descending =
    options.ags_groupset_order == AGSGroupsetOrder::DESCENDING or
    (options.ags_groupset_order == AGSGroupsetOrder::AUTO and options.adjoint);
  // The groupset solvers in the order of solution.
  auto ordered_solvers = wgs_solvers_;
  if (descending)
    std::reverse(ordered_solvers.begin(), ordered_solvers.end());
  const bool print_ags_stats = verbose_ and is_multi_groupset;
  last_solve_ = {};

  // The change between passes needs the scalar flux of the previous pass. It is kept only when it
  // is reported: with several groupsets, or when more than one pass is allowed. The source moments
  // need no copy, because each groupset solve leaves them as it found them.
  const bool measures_change = is_multi_groupset or max_iterations_ > 1;
  std::vector<double> phi_old;
  if (measures_change)
    phi_old = lbs_problem_.GetPhiOldLocal();

  // Without coupling to earlier groupsets, a second pass would only continue the groupset solves
  // that already converged, so the iteration stops after the first pass in which they all do.
  const auto inexact_reason =
    is_multi_groupset ? CheckAGSIterates(lbs_problem_, ordered_solvers) : std::nullopt;
  const bool one_pass_exact = is_multi_groupset and not inexact_reason;
  // A power iteration or time step solves many times with the same reason, so it is logged when
  // it changes.
  if (print_ags_stats and inexact_reason and *inexact_reason != logged_inexact_reason_)
    log.Log() << no_wrap << program_timer.GetTimeString() << " AGS iterates: " << *inexact_reason;
  logged_inexact_reason_ = inexact_reason.value_or("");

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
    for (auto& solver : ordered_solvers)
    {
      solver->Setup();
      solver->Solve();
      ValidateFiniteScalarFlux(lbs_problem_, solver);
      // Groupsets are solved one at a time, so only the groupset being solved needs its solver
      // storage. With several groupsets, every solve of a groupset, including those of each power
      // iteration and time step, allocates the storage in Setup() and frees it here. A single
      // groupset keeps its storage, since no other groupset needs the memory.
      if (is_multi_groupset)
        solver->ReleaseWorkspace();
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

    if (measures_change and lbs_problem_.GetOptions().ags_pointwise_convergence)
    {
      double pw_change = ComputePointwisePhiChange(lbs_problem_, phi_old);
      last_error = pw_change;
      const double rho = rate_estimator.Update(pw_change);

      AppendNumericField(iter_stats, "pw_change", pw_change, Scientific(6));
      AppendNumericField(iter_stats, "rho_est", rho, Fixed(4));

      if (rate_estimator.Converged(tolerance_))
        converged = true;
    }
    else if (measures_change)
    {
      double norm = ComputeL2PhiChange(lbs_problem_, phi_old);
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

    // A within-group solve that stops at its iteration limit is an inexact inner solve; the
    // across-groupset iteration continues and its own convergence test decides.
    if (wgs_status == IterationStatus::FAILED)
    {
      failed = true;
      break;
    }
    else if (converged)
      break;
    else if (measures_change)
      phi_old = lbs_problem_.GetPhiNewLocal();
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
