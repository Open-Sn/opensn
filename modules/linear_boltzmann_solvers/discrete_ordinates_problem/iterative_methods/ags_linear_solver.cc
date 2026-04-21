// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/iterative_methods/ags_linear_solver.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/iterative_methods/convergence.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/iterative_methods/sweep_wgs_context.h"
#include "modules/linear_boltzmann_solvers/lbs_problem/iterative_methods/iteration_logging.h"
#include "modules/linear_boltzmann_solvers/lbs_problem/lbs_problem.h"
#include "framework/mesh/mesh_continuum/mesh_continuum.h"
#include "framework/logging/log.h"
#include "framework/runtime.h"
#include "caliper/cali.h"

namespace opensn
{

void
AGSLinearSolver::Solve()
{
  CALI_CXX_MARK_SCOPE("AGSLinearSolver::Solve");
  const bool is_multi_groupset = wgs_solvers_.size() > 1;
  const bool print_ags_stats = verbose_ and is_multi_groupset;
  last_solve_ = {};

  phi_old_ = lbs_problem_.GetPhiOldLocal();

  // Save qmoms to be restored after each iteration. This is necessary for multiple ags iterations
  // to function and for keigen-value problems
  const auto saved_qmoms = lbs_problem_.GetQMomentsLocal();

  double pw_change_prev = 1.0;
  double last_error = 0.0;
  bool converged = false;
  bool failed = false;
  bool child_limited = false;
  unsigned int num_iterations = 0;
  const bool reports_as_iteration = is_multi_groupset;
  for (unsigned int iter = 0; iter < max_iterations_; ++iter)
  {
    num_iterations = iter + 1;
    for (auto& solver : wgs_solvers_)
    {
      solver->Setup();
      solver->Solve();
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
      double rho = (iter == 0) ? 0.0 : sqrt(pw_change / pw_change_prev);
      pw_change_prev = pw_change;

      AppendNumericField(iter_stats, "pw_change", pw_change, Scientific(6));
      AppendNumericField(iter_stats, "rho_est", rho, Fixed(4));

      if (pw_change < std::max(tolerance_ * (1.0 - rho), 1.0e-10))
        converged = true;
    }
    else
    {
      double norm = ComputeL2PhiChange(lbs_problem_, phi_old_);
      last_error = norm;

      AppendNumericField(iter_stats, "l2_change", norm, Scientific(6));

      if (norm < tolerance_)
        converged = true;
    }

    const auto iteration_status = IterationStatusFromSolve(converged, false, wgs_status);
    if (iteration_status != IterationStatus::NONE)
      iter_stats << ", status = " << IterationStatusName(iteration_status);
    iter_stats << FormatNestedStatusCounts("WGS", wgs_summaries);

    if (print_ags_stats)
      log.Log() << iter_stats.str();

    // Restore qmoms
    lbs_problem_.SetQMomentsFrom(saved_qmoms);

    if (wgs_status == IterationStatus::FAILED)
    {
      failed = true;
      break;
    }
    else if (wgs_status == IterationStatus::LIMIT)
    {
      child_limited = true;
      break;
    }
    else if (converged)
      break;
    else
      phi_old_ = lbs_problem_.GetPhiNewLocal();
  }

  const auto inner_status = failed
                              ? IterationStatus::FAILED
                              : (child_limited ? IterationStatus::LIMIT : IterationStatus::NONE);
  last_solve_.num_iterations = reports_as_iteration ? num_iterations : 0;
  last_solve_.status = reports_as_iteration
                         ? IterationStatusFromSolve(converged, true, inner_status)
                         : IterationStatus::NONE;
  last_solve_.metric_name =
    lbs_problem_.GetOptions().ags_pointwise_convergence ? "pw_change" : "l2_change";
  last_solve_.metric_value = last_error;
  if (print_ags_stats and not converged)
    log.Log() << program_timer.GetTimeString() << " " << FormatIterationSummary("AGS", last_solve_);

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
    const size_t num_angles = groupset.quadrature->abscissae.size();
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
    log.Log() << program_timer.GetTimeString() << " " << sweep_timing.str();

    std::stringstream sweep_work;
    sweep_work << label;
    AppendNumericField(sweep_work, "unknowns", num_unknowns, false);
    AppendNumericField(sweep_work, "lagged_unknowns", num_delayed_unknowns);
    AppendNumericField(sweep_work, "lagged_pct", delayed_unknown_percent, Fixed(2));
    log.Log() << program_timer.GetTimeString() << " " << sweep_work.str();

    sweep_context->ResetAccumulatedSweepStats();
  }
}

} // namespace opensn
