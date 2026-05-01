// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/iterative_methods/power_iteration_keigen.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/discrete_ordinates_problem.h"
#include "modules/linear_boltzmann_solvers/lbs_problem/lbs_problem.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/iterative_methods/ags_linear_solver.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/iterative_methods/wgs_context.h"
#include "modules/linear_boltzmann_solvers/lbs_problem/iterative_methods/iteration_logging.h"
#include "modules/linear_boltzmann_solvers/lbs_problem/compute/lbs_compute.h"
#include "framework/runtime.h"
#include "framework/logging/log.h"
#include "framework/utils/caliper_scopes.h"
#include "framework/utils/timer.h"
#include "caliper/cali.h"

namespace opensn
{

void
PowerIterationKEigenSolver(LBSProblem& lbs_problem,
                           double tolerance,
                           unsigned int max_iterations,
                           double& k_eff)
{
  CaliperPhaseScope cali_solve_phase("Solve", CaliperSolvePhaseDepth());
  CaliperRegionScope cali_pi("PI", CaliperPIScopeDepth());

  const std::string fname = "PowerIterationKEigenSolver";
  auto* do_problem = dynamic_cast<DiscreteOrdinatesProblem*>(&lbs_problem);
  if (not do_problem)
    throw std::logic_error(fname + ": requires a DiscreteOrdinatesProblem.");

  for (size_t gsid = 0; gsid < do_problem->GetNumWGSSolvers(); ++gsid)
  {
    auto wgs_solver = do_problem->GetWGSSolver(gsid);
    auto context = wgs_solver->GetContext();
    auto wgs_context = std::dynamic_pointer_cast<WGSContext>(context);

    if (not wgs_context)
      throw std::logic_error(fname + ": Cast failed.");

    wgs_context->lhs_src_scope = APPLY_WGS_SCATTER_SOURCES;
    wgs_context->rhs_src_scope = APPLY_AGS_SCATTER_SOURCES | APPLY_FIXED_SOURCES;
  }

  const auto& phi_old_local = lbs_problem.GetPhiOldLocal();
  const auto& phi_new_local = lbs_problem.GetPhiNewLocal();
  auto ags_solver = do_problem->GetAGSSolver();
  const auto& groupsets = lbs_problem.GetGroupsets();
  auto active_set_source_function = lbs_problem.GetActiveSetSourceFunction();

  k_eff = 1.0;
  double k_eff_prev = 1.0;
  double k_eff_change = 1.0;
  double F_prev = ComputeFissionProduction(lbs_problem, phi_old_local);

  // Start power iterations
  const bool print_ags_iters =
    lbs_problem.GetOptions().verbose_inner_iterations and groupsets.size() > 1;
  ags_solver->SetVerbosity(print_ags_iters);
  unsigned int nit = 0;
  bool converged = false;
  CALI_CXX_MARK_LOOP_BEGIN(pi_iteration, "PIIteration");
  while (nit < max_iterations)
  {
    CALI_CXX_MARK_LOOP_ITERATION(pi_iteration, nit);

    lbs_problem.ZeroQMoments();
    for (const auto& groupset : groupsets)
    {
      CALI_CXX_MARK_SCOPE("Source");
      active_set_source_function(groupset,
                                 lbs_problem.GetQMomentsLocal(),
                                 phi_old_local,
                                 APPLY_AGS_FISSION_SOURCES | APPLY_WGS_FISSION_SOURCES);
    }

    lbs_problem.ScaleQMoments(1.0 / k_eff);

    // This solves the inners for transport
    ags_solver->Solve();

    // Recompute k-eigenvalue
    double F_new = ComputeFissionProduction(lbs_problem, phi_new_local);
    k_eff = F_new / F_prev * k_eff;
    double reactivity = (k_eff - 1.0) / k_eff;

    // Check convergence, bookkeeping
    k_eff_change = fabs(k_eff - k_eff_prev) / k_eff;
    k_eff_prev = k_eff;
    F_prev = F_new;
    nit += 1;

    if (k_eff_change < std::max(tolerance, 1.0e-12))
      converged = true;

    // Print iteration summary
    if (lbs_problem.GetOptions().verbose_outer_iterations)
    {
      std::vector<IterationSummary> wgs_summaries;
      wgs_summaries.reserve(do_problem->GetNumWGSSolvers());
      for (size_t gsid = 0; gsid < do_problem->GetNumWGSSolvers(); ++gsid)
      {
        auto wgs_context =
          std::dynamic_pointer_cast<WGSContext>(do_problem->GetWGSSolver(gsid)->GetContext());
        if (wgs_context and HasIterationStatus(wgs_context->last_solve))
          wgs_summaries.emplace_back(wgs_context->last_solve);
      }

      const auto ags_summary = ags_solver->GetLastSolveSummary();
      const auto ags_status =
        HasIterationStatus(ags_summary) ? ags_summary.status : IterationStatus::NONE;
      const auto inner_status =
        MostSevereIterationStatus(ags_status, MostSevereIterationStatus(wgs_summaries));
      const auto outer_status =
        IterationStatusFromSolve(converged, nit >= max_iterations, inner_status);
      log.Log() << program_timer.GetTimeString() << " "
                << FormatKEigenOuterIteration("PI",
                                              nit,
                                              k_eff,
                                              k_eff_change,
                                              reactivity * 1.0e5,
                                              ags_summary,
                                              wgs_summaries,
                                              outer_status);
    }

    if (converged)
      break;
  } // for k iterations
  CALI_CXX_MARK_LOOP_END(pi_iteration);

  // Print summary
  size_t total_sweeps = 0;
  for (const auto& groupset : groupsets)
  {
    auto wgs_solver = do_problem->GetWGSSolver(groupset.id);
    auto wgs_context = std::dynamic_pointer_cast<WGSContext>(wgs_solver->GetContext());
    if (wgs_context)
      total_sweeps += wgs_context->counter_applications_of_inv_op;
  }

  log.Log() << program_timer.GetTimeString() << " "
            << FormatKEigenFinalSummary("PI",
                                        k_eff,
                                        k_eff_change,
                                        total_sweeps,
                                        "sweeps",
                                        IterationStatusFromSolve(converged, true));
}

} // namespace opensn
