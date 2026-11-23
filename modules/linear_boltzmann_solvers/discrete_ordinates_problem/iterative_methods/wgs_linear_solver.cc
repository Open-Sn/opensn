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
  aliases_flux_ = false;
}

void
WGSLinearSolver::Solve()
{
  CaliperPhaseScope cali_solve_phase("Solve", CaliperSolvePhaseDepth());
  CaliperRegionScope cali_wgs("WGS", CaliperWGSScopeDepth());
  CaliperRegionScope cali_method(WGSMethodName(method_), CaliperWGSMethodScopeDepth());
  // Vectors that view phi_old and phi_new are stale once either is reallocated.
  if (aliases_flux_)
  {
    const auto& do_problem = std::dynamic_pointer_cast<WGSContext>(context_ptr_)->do_problem;
    if (do_problem.GetPhiOldLocal().GetGeneration() != phi_old_generation_ or
        do_problem.GetPhiNewLocal().GetGeneration() != phi_new_generation_)
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

  // The solution and right-hand-side vectors are the groupset's flux block followed by its local
  // delayed angular unknowns. Without delayed unknowns, x_ and b_ view the groupset's owned blocks
  // of phi_old and phi_new. The last block can also hold the delayed unknowns in the tail after it,
  // which the problem reserves when it builds the solvers. The owned phi_new block is free while
  // the right-hand side is needed: the Krylov operator redirects phi_new to its output, and
  // phi_new is next written after the solve. The count of delayed unknowns is global and the block
  // order is the same on every rank, so every rank takes the same branch.
  auto gs_context_ptr = std::dynamic_pointer_cast<WGSContext>(context_ptr_);
  auto& groupset = gs_context_ptr->groupset;
  auto& do_problem = gs_context_ptr->do_problem;
  auto& phi_old = do_problem.GetPhiOldLocal();
  auto& phi_new = do_problem.GetPhiNewLocal();
  const auto block = do_problem.GetFluxBlock(groupset);
  const auto num_delayed_psi = groupset.angle_agg ? groupset.angle_agg->GetNumDelayedAngularDOFs()
                                                  : std::pair<std::size_t, std::size_t>{0, 0};
  aliases_flux_ = num_delayed_psi.second == 0 or block + 1 == phi_old.GetNumBlocks();
  const auto groupset_name = "WGSLinearSolver: groups [" + std::to_string(groupset.first_group) +
                             "-" + std::to_string(groupset.last_group) + "]";
  if (aliases_flux_)
  {
    const auto tail_size = num_delayed_psi.first;
    // The solver must not reallocate the flux storage, which others may be viewing.
    OpenSnLogicalErrorIf(tail_size > phi_old.GetTailCapacity() or
                           tail_size > phi_new.GetTailCapacity(),
                         groupset_name + " need " + std::to_string(tail_size) +
                           " tail entries for delayed angular unknowns on rank " +
                           std::to_string(opensn::mpi_comm.rank()) + ", but phi_old has " +
                           std::to_string(phi_old.GetTailCapacity()) + " and phi_new has " +
                           std::to_string(phi_new.GetTailCapacity()) + ".");
    const auto x_storage = phi_old.OwnedBlockWithTail(block, tail_size);
    const auto b_storage = phi_new.OwnedBlockWithTail(block, tail_size);
    OpenSnLogicalErrorIf(std::cmp_not_equal(x_storage.size(), num_local_dofs_),
                         groupset_name + " have " + std::to_string(x_storage.size()) +
                           " flux and delayed unknowns on rank " +
                           std::to_string(opensn::mpi_comm.rank()) + ", but the system has " +
                           std::to_string(num_local_dofs_) + ".");
    OpenSnPETScCall(VecCreateMPIWithArray(
      opensn::mpi_comm, 1, num_local_dofs_, num_global_dofs_, x_storage.data(), &x_));
    OpenSnPETScCall(VecCreateMPIWithArray(
      opensn::mpi_comm, 1, num_local_dofs_, num_global_dofs_, b_storage.data(), &b_));
    OpenSnPETScCall(VecSetOption(x_, VEC_IGNORE_NEGATIVE_INDICES, PETSC_TRUE));
    OpenSnPETScCall(VecSetOption(b_, VEC_IGNORE_NEGATIVE_INDICES, PETSC_TRUE));
    phi_old_generation_ = phi_old.GetGeneration();
    phi_new_generation_ = phi_new.GetGeneration();
  }
  else
  {
    x_ = CreateVector(num_local_dofs_, num_global_dofs_);
    OpenSnPETScCall(VecSet(x_, 0.0));
    OpenSnPETScCall(VecDuplicate(x_, &b_));
  }

  // The work vector is used only after the sources of the right-hand side are built, so it can use
  // the source scratch. Request its full size first, so that later redirects keep its address.
  const auto work_storage =
    gs_context_ptr->GetSourceScratch(static_cast<std::size_t>(num_local_dofs_));
  OpenSnPETScCall(VecCreateMPIWithArray(opensn::mpi_comm,
                                        1,
                                        num_local_dofs_,
                                        num_global_dofs_,
                                        work_storage.data(),
                                        &rhs_preconditioned_work_));
  OpenSnPETScCall(VecSetOption(rhs_preconditioned_work_, VEC_IGNORE_NEGATIVE_INDICES, PETSC_TRUE));

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

  if (aliases_flux_)
  {
    // The flux part of x_ is phi_old; append the delayed unknowns. Getting and restoring the array
    // also tells PETSc that x_ may have changed, since phi_old is written outside PETSc, so that
    // VecNorm does not return a norm it cached before.
    PetscScalar* x = nullptr;
    OpenSnPETScCall(VecGetArray(x_, &x));
    if (groupset.angle_agg)
    {
      const auto block_size =
        do_problem.GetPhiOldLocal().GetBlockSize(do_problem.GetFluxBlock(groupset));
      auto index = static_cast<int64_t>(block_size) - 1;
      groupset.angle_agg->AppendOldDelayedAngularDOFsToArray(index, x);
    }
    OpenSnPETScCall(VecRestoreArray(x_, &x));
  }
  else
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
  // When b_ views phi_new, this copy overwrites the right-hand side, which is no longer needed.
  const auto block = do_problem.GetFluxBlock(groupset);
  const auto block_size = do_problem.GetPhiOldLocal().GetBlockSize(block);
  const PetscScalar* x = nullptr;
  OpenSnPETScCall(VecGetArrayRead(x_, &x));
  if (not aliases_flux_)
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
