// SPDX-FileCopyrightText: 2025 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/iterative_methods/slepc_linear_keigen_solver.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/iterative_methods/ags_linear_solver.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/iterative_methods/wgs_context.h"
#include "modules/linear_boltzmann_solvers/lbs_problem/compute/lbs_compute.h"
#include "modules/linear_boltzmann_solvers/lbs_problem/vecops/lbs_vecops.h"
#include "framework/math/petsc_utils/petsc_utils.h"
#include "framework/math/math.h"
#include "framework/logging/log.h"
#include "framework/utils/timer.h"
#include "framework/runtime.h"
#include <algorithm>
#include <cassert>
#include <cmath>
#include <exception>
#include <iomanip>
#include <sstream>
#include <slepceps.h>

namespace opensn
{

namespace
{

struct SavedWGSSourceScopes
{
  std::vector<WGSContext*> contexts;
  std::vector<SourceFlags> lhs_scopes;
  std::vector<SourceFlags> rhs_scopes;

  explicit SavedWGSSourceScopes(DiscreteOrdinatesProblem& do_problem)
  {
    const auto num_wgs = do_problem.GetNumWGSSolvers();
    contexts.reserve(num_wgs);
    lhs_scopes.reserve(num_wgs);
    rhs_scopes.reserve(num_wgs);

    for (size_t i = 0; i < num_wgs; ++i)
    {
      auto wgs_solver = do_problem.GetWGSSolver(i);
      OpenSnLogicalErrorIf(not wgs_solver,
                           do_problem.GetName() + ": Null WGS solver in SLEPc k-eigen solve.");

      auto context = wgs_solver->GetContext();
      auto wgs_context = std::dynamic_pointer_cast<WGSContext>(context);
      OpenSnLogicalErrorIf(not wgs_context, do_problem.GetName() + ": Cast to WGSContext failed.");

      contexts.push_back(wgs_context.get());
      lhs_scopes.push_back(wgs_context->lhs_src_scope);
      rhs_scopes.push_back(wgs_context->rhs_src_scope);

      wgs_context->lhs_src_scope.Unset(APPLY_WGS_FISSION_SOURCES);
      wgs_context->rhs_src_scope.Unset(APPLY_AGS_FISSION_SOURCES);
    }
  }

  ~SavedWGSSourceScopes()
  {
    for (size_t i = 0; i < contexts.size(); ++i)
    {
      contexts[i]->lhs_src_scope = lhs_scopes[i];
      contexts[i]->rhs_src_scope = rhs_scopes[i];
    }
  }
};

struct SavedTransportState
{
  struct AngleAggregationState
  {
    AngleAggregation* angle_agg = nullptr;
    std::vector<double> psi_old;
    std::vector<double> psi_new;
  };

  DiscreteOrdinatesProblem& do_problem;
  std::vector<double> phi_old;
  std::vector<double> phi_new;
  std::vector<double> q_moments;
  std::vector<AngleAggregationState> angle_agg_states;

  explicit SavedTransportState(DiscreteOrdinatesProblem& do_problem)
    : do_problem(do_problem),
      phi_old(do_problem.GetPhiOldLocal()),
      phi_new(do_problem.GetPhiNewLocal()),
      q_moments(do_problem.GetQMomentsLocal())
  {
    angle_agg_states.reserve(do_problem.GetGroupsets().size());
    for (const auto& groupset : do_problem.GetGroupsets())
    {
      if (not groupset.angle_agg)
        continue;
      angle_agg_states.push_back({groupset.angle_agg.get(),
                                  groupset.angle_agg->GetOldDelayedAngularDOFsAsSTLVector(),
                                  groupset.angle_agg->GetNewDelayedAngularDOFsAsSTLVector()});
    }
  }

  ~SavedTransportState() noexcept
  {
    try
    {
      Restore();
    }
    catch (...)
    {
      std::terminate();
    }
  }

  void Restore()
  {
    CopyInto(do_problem.GetPhiOldLocal(), phi_old);
    CopyInto(do_problem.GetPhiNewLocal(), phi_new);
    CopyInto(do_problem.GetQMomentsLocal(), q_moments);
    for (const auto& state : angle_agg_states)
    {
      state.angle_agg->SetOldDelayedAngularDOFsFromSTLVector(state.psi_old);
      state.angle_agg->SetNewDelayedAngularDOFsFromSTLVector(state.psi_new);
    }
  }

  static void CopyInto(std::vector<double>& destination, const std::vector<double>& source)
  {
    assert(destination.size() == source.size() && "SavedTransportState restore size mismatch.");
    std::copy(source.begin(), source.end(), destination.begin());
  }
};

PetscErrorCode
ApplyLinearKEigenOperator(SLEPcLinearKEigenContext& ctx, Vec x, Vec y)
{
  auto& do_problem = ctx.do_problem;
  auto& phi_old_local = do_problem->GetPhiOldLocal();
  auto& phi_new_local = do_problem->GetPhiNewLocal();
  auto& q_moments_local = do_problem->GetQMomentsLocal();
  auto active_set_source_function = do_problem->GetActiveSetSourceFunction();
  auto ags_solver = do_problem->GetAGSSolver();

  LBSVecOps::SetPrimarySTLvectorFromGroupScopedPETScVec(
    *do_problem, 0, do_problem->GetNumGroups() - 1, x, PhiSTLOption::PHI_OLD);

  Set(q_moments_local, 0.0);

  // Apply the linear k-eigenvalue operator:
  //   A phi = (I - D L^-1 M S)^-1 D L^-1 F phi.
  // The fission source is the explicit input to the operator. The AGS solve applies the
  // scattering inverse with fission removed from the WGS source scopes.
  for (const auto& groupset : do_problem->GetGroupsets())
    active_set_source_function(groupset,
                               q_moments_local,
                               phi_old_local,
                               APPLY_AGS_FISSION_SOURCES | APPLY_WGS_FISSION_SOURCES);

  OpenSnLogicalErrorIf(not ags_solver, "SLEPcLinearKEigenSolver: AGS solver is not available.");
  ags_solver->Solve();

  LBSVecOps::SetGroupScopedPETScVecFromPrimarySTLvector(
    *do_problem, 0, do_problem->GetNumGroups() - 1, y, phi_new_local);

  return 0;
}

// Shell matrix multiplication: y = A * x
PetscErrorCode
ShellMult(Mat M, Vec x, Vec y)
{
  void* raw_ctx = nullptr;
  MatShellGetContext(M, static_cast<void*>(&raw_ctx));

  auto* ctx = static_cast<SLEPcLinearKEigenContext*>(raw_ctx);
  SavedTransportState saved_transport_state(*ctx->do_problem);

  return ApplyLinearKEigenOperator(*ctx, x, y);
}

} // namespace

static PetscErrorCode
SLEPcLinearKEigenMonitor(EPS unused_eps,
                         PetscInt its,
                         PetscInt unused_nconv,
                         PetscScalar* unused_eigr, // NOLINT(readability-non-const-parameter)
                         PetscScalar* unused_eigi, // NOLINT(readability-non-const-parameter)
                         PetscReal* errest,
                         PetscInt nest,
                         void* ctx)
{
  (void)unused_eps;
  (void)unused_nconv;
  (void)unused_eigr;
  (void)unused_eigi;

  const bool has_error_estimate = nest > 0 and errest != nullptr;
  if (ctx)
  {
    auto* kctx = static_cast<SLEPcLinearKEigenContext*>(ctx);
    if (has_error_estimate)
      kctx->last_eps_residual = errest[0];
  }
  std::stringstream iter_info;
  iter_info << program_timer.GetTimeString() << " SLEPc EPS iteration = " << its;
  if (has_error_estimate)
    iter_info << ", residual = " << std::scientific << std::setprecision(6) << errest[0]
              << std::defaultfloat;
  else
    iter_info << ", status = iterating";
  log.Log() << iter_info.str();

  return 0;
}

void
SLEPcLinearKEigenSolver::SetMonitor()
{
}

void
SLEPcLinearKEigenSolver::SetSystemSize()
{
  auto ctx = std::dynamic_pointer_cast<SLEPcLinearKEigenContext>(context_ptr_);
  // SLEPc solves the fission-production eigenproblem in scalar-flux space. Delayed angular
  // fluxes are still solved as part of each transport inverse application, but they are internal
  // operator state rather than EPS eigenvector entries.
  auto sizes = ctx->do_problem->LBSProblem::GetNumPhiIterativeUnknowns();
  num_local_dofs_ = static_cast<int64_t>(sizes.first);
  num_global_dofs_ = static_cast<int64_t>(sizes.second);
}

void
SLEPcLinearKEigenSolver::SetSystem()
{
  // Create shell matrix A
  MatCreateShell(PETSC_COMM_WORLD,
                 static_cast<PetscInt>(num_local_dofs_),
                 static_cast<PetscInt>(num_local_dofs_),
                 static_cast<PetscInt>(num_global_dofs_),
                 static_cast<PetscInt>(num_global_dofs_),
                 context_ptr_.get(),
                 &A_);
  // PETSc erases the shell-operation type to PetscErrorCodeFn*, so a cast is required here.
  // NOLINTNEXTLINE(cppcoreguidelines-pro-type-cstyle-cast,modernize-redundant-void-arg)
  MatShellSetOperation(A_, MATOP_MULT, (void (*)(void))ShellMult);

  // Create x
  x_ = CreateVector(num_local_dofs_, num_global_dofs_);
}

void
SLEPcLinearKEigenSolver::SetInitialGuess()
{
  auto ctx = std::dynamic_pointer_cast<SLEPcLinearKEigenContext>(context_ptr_);
  LBSVecOps::SetGroupScopedPETScVecFromPrimarySTLvector(*ctx->do_problem,
                                                        0,
                                                        ctx->do_problem->GetNumGroups() - 1,
                                                        x_,
                                                        ctx->do_problem->GetPhiOldLocal());
}

void
SLEPcLinearKEigenSolver::Solve()
{
  PreSolveCallback();
  auto ctx = std::dynamic_pointer_cast<SLEPcLinearKEigenContext>(context_ptr_);
  SavedWGSSourceScopes saved_wgs_scopes(*ctx->do_problem);
  SetInitialGuess();

  EPSSetFromOptions(eps_);
  EPSSetOperators(eps_, A_, nullptr);
  EPSSetProblemType(eps_, EPS_NHEP);
  EPSSetWhichEigenpairs(eps_, EPS_LARGEST_MAGNITUDE);
  EPSSetTolerances(eps_, tolerance_options.residual_absolute, tolerance_options.maximum_iterations);
  EPSSetType(eps_, eps_type_.c_str());
  EPSSetInitialSpace(eps_, 1, &x_);
  EPSMonitorSet(eps_, SLEPcLinearKEigenMonitor, context_ptr_.get(), nullptr);
  EPSSolve(eps_);

  // Check for convergence
  PetscInt nconv = 0;
  EPSGetConverged(eps_, &nconv);
  if (nconv == 0)
  {
    OpenSnLogicalError(program_timer.GetTimeString() + " SLEPc EPS: No eigenpairs converged.");
  }

  PostSolveCallback();
}

void
SLEPcLinearKEigenSolver::PreSolveCallback()
{
  log.Log() << "Executing SLEPc Eigenvalue Problem Solver with iterative method "
            << GetIterativeMethodName() << std::endl;

  auto ctx = std::dynamic_pointer_cast<SLEPcLinearKEigenContext>(context_ptr_);
  ctx->last_eps_residual = 1.0;
  ctx->last_operator_residual = 1.0;
}

void
SLEPcLinearKEigenSolver::PostSolveCallback()
{
  auto ctx = std::dynamic_pointer_cast<SLEPcLinearKEigenContext>(context_ptr_);
  auto& do_problem = ctx->do_problem;
  auto& phi_old_local = do_problem->GetPhiOldLocal();
  auto& phi_new_local = do_problem->GetPhiNewLocal();

  // Extract the converged eigenvector and compute k as the fission-production ratio.
  PetscScalar kr = 0.0;
  PetscScalar ki = 0.0;
  Vec x = CreateVector(num_local_dofs_, num_global_dofs_);
  EPSGetEigenpair(eps_, 0, &kr, &ki, x, nullptr);
  LBSVecOps::SetPrimarySTLvectorFromGroupScopedPETScVec(
    *do_problem, 0, do_problem->GetNumGroups() - 1, x, PhiSTLOption::PHI_OLD);
  double F_prev = ComputeFissionProduction(*do_problem, phi_old_local);

  // Perform a transport solve, get the resulting eigenvector as phi_new, and compute fission rate
  Vec y = CreateVector(num_local_dofs_, num_global_dofs_);
  ApplyLinearKEigenOperator(*ctx, x, y);
  LBSVecOps::SetPrimarySTLvectorFromGroupScopedPETScVec(
    *do_problem, 0, do_problem->GetNumGroups() - 1, y, PhiSTLOption::PHI_NEW);
  double F = ComputeFissionProduction(*do_problem, phi_new_local);

  // Compute k-eigenvalue
  OpenSnLogicalErrorIf(F_prev == 0.0,
                       "SLEPcLinearKEigenSolver: previous fission production is zero. "
                       "Cannot compute k-eigenvalue.");
  ctx->eigenvalue = F / F_prev;
  log.Log() << "Final k-eigenvalue: " << std::setprecision(7) << ctx->eigenvalue << "\n";

  // Compute operator residual ||A x - k x|| / ||x|| using the same operator SLEPc solved.
  Vec residual = CreateVector(num_local_dofs_, num_global_dofs_);
  OpenSnPETScCall(VecCopy(y, residual));
  OpenSnPETScCall(VecAXPY(residual, -ctx->eigenvalue, x));
  PetscReal x_norm = 0.0;
  PetscReal residual_norm = 0.0;
  OpenSnPETScCall(VecNorm(x, NORM_2, &x_norm));
  OpenSnPETScCall(VecNorm(residual, NORM_2, &residual_norm));
  ctx->last_operator_residual = (x_norm > 0.0) ? residual_norm / x_norm : residual_norm;
  log.Log() << "Operator residual: " << std::setprecision(7) << ctx->last_operator_residual << "\n";

  OpenSnLogicalErrorIf(not std::isfinite(ctx->last_operator_residual),
                       "SLEPcLinearKEigenSolver: final operator residual is not finite.");

  const auto ags_solver = do_problem->GetAGSSolver();
  const double residual_tolerance = std::max({10.0 * tolerance_options.residual_absolute,
                                              ags_solver ? 10.0 * ags_solver->GetTolerance() : 0.0,
                                              1.0e-10});
  if (ctx->last_operator_residual > residual_tolerance)
  {
    std::stringstream warning;
    warning << "SLEPcLinearKEigenSolver: final operator residual " << std::scientific
            << std::setprecision(6) << ctx->last_operator_residual
            << " exceeds diagnostic tolerance " << residual_tolerance
            << ". This residual includes the accuracy of the inner transport solve.";
    log.Log0Warning() << warning.str();
  }

  // Leave the problem state consistent with the operator application.
  LBSVecOps::SetPrimarySTLvectorFromGroupScopedPETScVec(
    *do_problem, 0, do_problem->GetNumGroups() - 1, x, PhiSTLOption::PHI_OLD);
  LBSVecOps::SetPrimarySTLvectorFromGroupScopedPETScVec(
    *do_problem, 0, do_problem->GetNumGroups() - 1, y, PhiSTLOption::PHI_NEW);

  VecDestroy(&x);
  VecDestroy(&y);
  VecDestroy(&residual);
}

} // namespace opensn
