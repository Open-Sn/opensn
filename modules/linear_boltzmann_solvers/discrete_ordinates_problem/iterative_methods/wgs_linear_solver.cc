// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/iterative_methods/wgs_linear_solver.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/iterative_methods/wgs_convergence_test.h"
#include "modules/linear_boltzmann_solvers/lbs_problem/vecops/lbs_vecops.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/discrete_ordinates_problem.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/iterative_methods/sweep_wgs_context.h"
#include "framework/math/petsc_utils/petsc_utils.h"
#include "framework/logging/log.h"
#include "framework/utils/caliper_scopes.h"
#include "framework/utils/error.h"
#include "framework/utils/timer.h"
#include "framework/runtime.h"
#include <petscksp.h>
#include "caliper/cali.h"
#include <algorithm>
#include <memory>
#include <span>
#include <string>
#include <utility>

namespace opensn
{

namespace
{

const char*
WGSMethodName(LinearSystemSolver::IterativeMethod method)
{
  switch (method)
  {
    case LinearSystemSolver::IterativeMethod::PETSC_RICHARDSON:
      return "Richardson";
    case LinearSystemSolver::IterativeMethod::PETSC_GMRES:
      return "GMRES";
    case LinearSystemSolver::IterativeMethod::PETSC_FGMRES:
      return "FGMRES";
    case LinearSystemSolver::IterativeMethod::PETSC_BICGSTAB:
      return "BiCGSTAB";
    default:
      return "WGSIteration";
  }
}

} // namespace

WGSLinearSolver::WGSLinearSolver(const std::shared_ptr<WGSContext>& gs_context_ptr)
  : PETScLinearSolver(gs_context_ptr->groupset.iterative_method, gs_context_ptr),
    rhs_preconditioned_work_(nullptr)
{
  auto& groupset = gs_context_ptr->groupset;
  auto& solver_tol_options = this->GetToleranceOptions();
  solver_tol_options.residual_absolute = groupset.residual_tolerance;
  solver_tol_options.maximum_iterations = static_cast<PetscInt>(groupset.max_iterations);
  solver_tol_options.gmres_restart_interval = static_cast<PetscInt>(groupset.gmres_restart_intvl);
}

WGSLinearSolver::~WGSLinearSolver()
{
  OpenSnPETScCall(VecDestroy(&rhs_preconditioned_work_));
}

void
WGSLinearSolver::ReleaseWorkspace()
{
  OpenSnPETScCall(VecDestroy(&rhs_preconditioned_work_));
  std::dynamic_pointer_cast<WGSContext>(context_ptr_)->ReleaseSourceScratch();
  PETScLinearSolver::ReleaseWorkspace();
  x_aliases_phi_old_ = false;
}

void
WGSLinearSolver::Solve()
{
  CaliperPhaseScope cali_solve_phase("Solve", CaliperSolvePhaseDepth());
  CaliperRegionScope cali_wgs("WGS", CaliperWGSScopeDepth());
  CaliperRegionScope cali_method(WGSMethodName(method_), CaliperWGSMethodScopeDepth());
  // A solution vector that views phi_old is stale once phi_old is reallocated.
  if (x_aliases_phi_old_)
  {
    const auto& do_problem = std::dynamic_pointer_cast<WGSContext>(context_ptr_)->do_problem;
    if (do_problem.GetPhiOldLocal().GetGeneration() != x_alias_generation_)
      ReleaseWorkspace();
  }
  PETScLinearSolver::Solve();
}

void
WGSLinearSolver::PreSetupCallback()
{
  auto gs_context_ptr = std::dynamic_pointer_cast<WGSContext>(context_ptr_);
  gs_context_ptr->PreSetupCallback();
}

void
WGSLinearSolver::SetConvergenceTest()
{
  OpenSnPETScCall(KSPSetConvergenceTest(ksp_, &GSConvergenceTest, nullptr, nullptr));
}

void
WGSLinearSolver::SetSystemSize()
{
  auto gs_context_ptr = std::dynamic_pointer_cast<WGSContext>(context_ptr_);
  const auto sizes = gs_context_ptr->GetSystemSize();
  num_local_dofs_ = sizes.first;
  num_global_dofs_ = sizes.second;
}

void
WGSLinearSolver::SetSystem()
{
  if (IsSystemSet())
    return;

  // Without delayed angular unknowns, the solution vector is exactly the groupset's phi_old block,
  // so x_ views that block and the initial guess and the solution need no copies. The count of
  // delayed unknowns is global, so every rank takes the same branch.
  auto gs_context_ptr = std::dynamic_pointer_cast<WGSContext>(context_ptr_);
  auto& groupset = gs_context_ptr->groupset;
  auto& phi_old = gs_context_ptr->do_problem.GetPhiOldLocal();
  x_aliases_phi_old_ =
    not groupset.angle_agg or groupset.angle_agg->GetNumDelayedAngularDOFs().second == 0;
  if (x_aliases_phi_old_)
  {
    const auto owned = phi_old.OwnedBlock(gs_context_ptr->do_problem.GetFluxBlock(groupset));
    OpenSnLogicalErrorIf(
      std::cmp_not_equal(owned.size(), num_local_dofs_),
      "WGSLinearSolver: the flux block of groups [" + std::to_string(groupset.first_group) + "-" +
        std::to_string(groupset.last_group) + "] has " + std::to_string(owned.size()) +
        " entries on rank " + std::to_string(opensn::mpi_comm.rank()) + ", but the system has " +
        std::to_string(num_local_dofs_) + ".");
    OpenSnPETScCall(VecCreateMPIWithArray(
      opensn::mpi_comm, 1, num_local_dofs_, num_global_dofs_, owned.data(), &x_));
    OpenSnPETScCall(VecSetOption(x_, VEC_IGNORE_NEGATIVE_INDICES, PETSC_TRUE));
    x_alias_generation_ = phi_old.GetGeneration();
  }
  else
  {
    x_ = CreateVector(num_local_dofs_, num_global_dofs_);
    OpenSnPETScCall(VecSet(x_, 0.0));
  }
  OpenSnPETScCall(VecDuplicate(x_, &b_));
  OpenSnPETScCall(VecDuplicate(x_, &rhs_preconditioned_work_));

  // Create the matrix-shell
  OpenSnPETScCall(MatCreateShell(opensn::mpi_comm,
                                 num_local_dofs_,
                                 num_local_dofs_,
                                 num_global_dofs_,
                                 num_global_dofs_,
                                 &(*context_ptr_),
                                 &A_));

  // Set the action-operator
  // NOLINTNEXTLINE(cppcoreguidelines-pro-type-cstyle-cast,modernize-redundant-void-arg)
  OpenSnPETScCall(MatShellSetOperation(A_, MATOP_MULT, (void (*)(void))LinearSolverMatrixAction));

  // Set solver operators
  OpenSnPETScCall(KSPSetOperators(ksp_, A_, A_));
  OpenSnPETScCall(KSPSetUp(ksp_));
}

void
WGSLinearSolver::SetPreconditioner()
{
  if (IsSystemSet())
    return;
  auto gs_context_ptr = std::dynamic_pointer_cast<WGSContext>(context_ptr_);
  gs_context_ptr->SetPreconditioner(ksp_);
}

void
WGSLinearSolver::PostSetupCallback()
{
  auto gs_context_ptr = std::dynamic_pointer_cast<WGSContext>(context_ptr_);
  gs_context_ptr->PostSetupCallback();
}

void
WGSLinearSolver::PreSolveCallback()
{
  auto gs_context_ptr = std::dynamic_pointer_cast<WGSContext>(context_ptr_);
  gs_context_ptr->last_solve = {};
  gs_context_ptr->PreSolveCallback();
}

void
WGSLinearSolver::SetInitialGuess()
{
  // The current flux of the groupset is the initial guess unless it is zero, in which case PETSc
  // can skip the operator application that a nonzero guess needs. The flux has the scale of the
  // source, so any nonzero norm is meaningful: a solution whose norm happens to be small is still
  // the best guess, and discarding it restarts the solve from zero.

  auto gs_context_ptr = std::dynamic_pointer_cast<WGSContext>(context_ptr_);
  auto& groupset = gs_context_ptr->groupset;
  auto& do_problem = gs_context_ptr->do_problem;

  if (not x_aliases_phi_old_)
    LBSVecOps::SetGSPETScVecFromPrimarySTLvector(do_problem, groupset, x_, PhiSTLOption::PHI_OLD);

  double init_guess_norm = 0.0;
  OpenSnPETScCall(VecNorm(x_, NORM_2, &init_guess_norm));
  OpenSnPETScCall(KSPSetInitialGuessNonzero(ksp_, PETSC_FALSE));

  if (init_guess_norm > 0.0)
    OpenSnPETScCall(KSPSetInitialGuessNonzero(ksp_, PETSC_TRUE));
}

void
WGSLinearSolver::SetRHS()
{
  CALI_CXX_MARK_SCOPE("SetRHS");

  auto gs_context_ptr = std::dynamic_pointer_cast<WGSContext>(context_ptr_);
  auto& groupset = gs_context_ptr->groupset;
  auto& do_problem = gs_context_ptr->do_problem;

  // If we have a single richardson iteration then the user probably wants only a single sweep.
  // Therefore, we are going to combine the scattering source (normally included in the
  // lhs_src_scope) into the sweep for the RHS, sweep into the solution, and just suppress the
  // kspsolve part.
  const bool single_richardson =
    groupset.iterative_method == LinearSystemSolver::IterativeMethod::PETSC_RICHARDSON and
    tolerance_options.maximum_iterations == 1;
  const auto scope = single_richardson
                       ? gs_context_ptr->GetRHSSourceScope() | gs_context_ptr->GetLHSSourceScope()
                       : gs_context_ptr->GetRHSSourceScope() | ZERO_INCOMING_DELAYED_PSI;
  Vec rhs = single_richardson ? x_ : b_;

  auto& phi_new = do_problem.GetPhiNewLocal();
  const auto block = do_problem.GetFluxBlock(groupset);
  const auto block_size = phi_new.GetBlockSize(block);
  PetscScalar* rhs_array = nullptr;
  OpenSnPETScCall(VecGetArray(rhs, &rhs_array));
  {
    // Sources go to scratch, so q is left unchanged, and the sweep writes the flux part of the
    // right-hand side directly.
    const auto q_override = gs_context_ptr->RedirectSourceMoments(false);
    {
      CALI_CXX_MARK_SCOPE("Source");
      gs_context_ptr->set_source_function(groupset,
                                          do_problem.GetQMomentsLocal().Blocks(),
                                          do_problem.GetPhiOldLocal().ConstBlocks(),
                                          scope);
    }

    // Enable RHS time (tau*psi^n)
    auto sweep_ctx = std::dynamic_pointer_cast<SweepWGSContext>(gs_context_ptr);
    if (not single_richardson and sweep_ctx and sweep_ctx->sweep_chunk->IsTimeDependent())
      sweep_ctx->sweep_chunk->IncludeRHSTimeTerm(true);

    // Apply transport operator
    const GroupsetBlockedVector::ScopedBlockOverride phi_new_override(
      phi_new, block, std::span<double>(rhs_array, block_size));
    gs_context_ptr->ApplyInverseTransportOperator(scope);

    if (groupset.angle_agg)
    {
      auto index = static_cast<int64_t>(block_size) - 1;
      groupset.angle_agg->AppendNewDelayedAngularDOFsToArray(index, rhs_array);
    }
  }
  OpenSnPETScCall(VecRestoreArray(rhs, &rhs_array));

  // Compute RHS norm
  OpenSnPETScCall(VecNorm(rhs, NORM_2, &gs_context_ptr->rhs_norm));

  // Compute precondition RHS norm
  PC pc = nullptr;
  OpenSnPETScCall(KSPGetPC(ksp_, &pc));
  OpenSnPETScCall(PCApply(pc, rhs, rhs_preconditioned_work_));
  OpenSnPETScCall(
    VecNorm(rhs_preconditioned_work_, NORM_2, &gs_context_ptr->rhs_preconditioned_norm));

  if (single_richardson)
    SetKSPSolveSuppressionFlag(true);
}

void
WGSLinearSolver::PostSolveCallback()
{
  auto gs_context_ptr = std::dynamic_pointer_cast<WGSContext>(context_ptr_);
  auto& groupset = gs_context_ptr->groupset;
  auto& do_problem = gs_context_ptr->do_problem;
  const bool report_nonfatal = do_problem.GetOptions().verbose_inner_iterations;

  // Get convergence reason
  if (not GetKSPSolveSuppressionFlag())
  {
    KSPConvergedReason reason = KSP_CONVERGED_ITERATING;
    OpenSnPETScCall(KSPGetConvergedReason(ksp_, &reason));
    PetscInt its = 0;
    OpenSnPETScCall(KSPGetIterationNumber(ksp_, &its));
    gs_context_ptr->last_solve.num_iterations = static_cast<unsigned int>(its);
    gs_context_ptr->last_solve.status = KSPReasonToIterationStatus(reason);
    gs_context_ptr->last_solve.detail = GetPETScConvergedReasonstring(reason);
    if (report_nonfatal and (gs_context_ptr->last_solve.status == IterationStatus::FAILED or
                             gs_context_ptr->last_solve.status == IterationStatus::LIMIT))
    {
      log.Log() << no_wrap << program_timer.GetTimeString() << " "
                << FormatIterationSummary("WGS groups [" + std::to_string(groupset.first_group) +
                                            "-" + std::to_string(groupset.last_group) + "]",
                                          gs_context_ptr->last_solve);
    }
  }

  // Copy x to the local solution. The flux part of x is the groupset's block in order. A view of
  // phi_old already holds it, and a context that rebuilds phi_new after the solve overwrites it.
  const auto block = do_problem.GetFluxBlock(groupset);
  const auto block_size = do_problem.GetPhiOldLocal().GetBlockSize(block);
  const PetscScalar* x = nullptr;
  OpenSnPETScCall(VecGetArrayRead(x_, &x));
  if (not x_aliases_phi_old_)
    std::copy_n(x, block_size, do_problem.GetPhiOldLocal().Block(block).begin());
  if (not gs_context_ptr->RebuildsFluxAfterSolve())
    std::copy_n(x, block_size, do_problem.GetPhiNewLocal().Block(block).begin());
  if (groupset.angle_agg)
  {
    auto index = static_cast<int64_t>(block_size) - 1;
    groupset.angle_agg->SetNewDelayedAngularDOFsFromArray(index, x);
    index = static_cast<int64_t>(block_size) - 1;
    groupset.angle_agg->SetOldDelayedAngularDOFsFromArray(index, x);
  }
  OpenSnPETScCall(VecRestoreArrayRead(x_, &x));

  // Context specific callback
  gs_context_ptr->PostSolveCallback();
}

} // namespace opensn
