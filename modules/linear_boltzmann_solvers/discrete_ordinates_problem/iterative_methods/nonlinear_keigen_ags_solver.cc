// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/iterative_methods/nonlinear_keigen_ags_solver.h"
#include "modules/linear_boltzmann_solvers/lbs_problem/iterative_methods/iteration_logging.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/iterative_methods/snes_k_monitor.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/iterative_methods/nonlinear_keigen_ags_residual_func.h"
#include "modules/linear_boltzmann_solvers/lbs_problem/vecops/lbs_vecops.h"
#include "modules/linear_boltzmann_solvers/lbs_problem/compute/lbs_compute.h"
#include "framework/math/petsc_utils/petsc_utils.h"
#include "framework/runtime.h"
#include "framework/logging/log.h"
#include "framework/utils/timer.h"
#include <petscsnes.h>

namespace opensn
{

namespace
{

std::shared_ptr<NLKEigenAGSContext>
GetNLKAGSContextPtr(const std::shared_ptr<NonLinearSolverContext>& context,
                    const std::string& func_name)
{
  auto nlk_ags_context = std::dynamic_pointer_cast<NLKEigenAGSContext>(context);
  if (not nlk_ags_context)
    throw std::runtime_error(func_name + ": context casting failure");
  return nlk_ags_context;
}

} // namespace

void
NLKEigenvalueAGSSolver::PreSetupCallback()
{
  auto nl_context_ptr = GetNLKAGSContextPtr(context_ptr_, __PRETTY_FUNCTION__);

  auto& lbs_problem = nl_context_ptr->lbs_problem;
  nl_context_ptr->groupset_ids.clear();
  for (const auto& groupset : lbs_problem->GetGroupsets())
  {
    if ((groupset.apply_wgdsa or groupset.apply_tgdsa) and lbs_problem->GetNumGroupsets() > 1)
      throw std::logic_error(
        std::string(__PRETTY_FUNCTION__) +
        ": Preconditioning currently only supports single groupset simulations.");
    nl_context_ptr->groupset_ids.push_back(groupset.id);
  }
}

void
NLKEigenvalueAGSSolver::SetMonitor()
{
  auto nl_context_ptr = GetNLKAGSContextPtr(context_ptr_, __PRETTY_FUNCTION__);

  auto& lbs_problem = nl_context_ptr->lbs_problem;
  if (lbs_problem->GetOptions().verbose_outer_iterations)
    OpenSnPETScCall(SNESMonitorSet(
      nl_solver_, &KEigenSNESMonitor, &nl_context_ptr->kresid_func_context, nullptr));

  if (lbs_problem->GetOptions().verbose_inner_iterations)
  {
    KSP ksp = nullptr;
    OpenSnPETScCall(SNESGetKSP(nl_solver_, &ksp));
    OpenSnPETScCall(
      KSPMonitorSet(ksp, &KEigenKSPMonitor, &nl_context_ptr->kresid_func_context, nullptr));
  }
}

void
NLKEigenvalueAGSSolver::SetSystemSize()
{
  auto nl_context_ptr = GetNLKAGSContextPtr(context_ptr_, __PRETTY_FUNCTION__);

  auto& lbs_problem = nl_context_ptr->lbs_problem;
  auto sizes = lbs_problem->GetNumPhiIterativeUnknowns();

  num_local_dofs_ = static_cast<PetscInt>(sizes.first);
  num_global_dofs_ = static_cast<PetscInt>(sizes.second);
}

void
NLKEigenvalueAGSSolver::SetSystem()
{
  // Create the vectors
  x_ = CreateVector(num_local_dofs_, num_global_dofs_);
  OpenSnPETScCall(VecDuplicate(x_, &r_));
}

void
NLKEigenvalueAGSSolver::SetFunction()
{
  auto nl_context_ptr = GetNLKAGSContextPtr(context_ptr_, __PRETTY_FUNCTION__);

  OpenSnPETScCall(SNESSetFunction(
    nl_solver_, r_, NLKEigenResidualFunction, &nl_context_ptr->kresid_func_context));
}

void
NLKEigenvalueAGSSolver::SetJacobian()
{
  OpenSnPETScCall(MatCreateSNESMF(nl_solver_, &J_));
  OpenSnPETScCall(SNESSetJacobian(nl_solver_, J_, J_, MatMFFDComputeJacobian, nullptr));
}

void
NLKEigenvalueAGSSolver::SetInitialGuess()
{
  auto nl_context_ptr = GetNLKAGSContextPtr(context_ptr_, __PRETTY_FUNCTION__);

  auto& lbs_problem = nl_context_ptr->lbs_problem;
  const auto& groupset_ids = nl_context_ptr->groupset_ids;

  LBSVecOps::SetMultiGSPETScVecFromPrimarySTLvector(
    *lbs_problem, groupset_ids, x_, PhiSTLOption::PHI_OLD);
}

void
NLKEigenvalueAGSSolver::PostSolveCallback()
{
  auto nl_context_ptr = GetNLKAGSContextPtr(context_ptr_, __PRETTY_FUNCTION__);

  auto& lbs_problem = nl_context_ptr->lbs_problem;
  const auto& groupset_ids = nl_context_ptr->groupset_ids;

  // Unpack solution
  LBSVecOps::SetPrimarySTLvectorFromMultiGSPETScVec(
    *lbs_problem, groupset_ids, x_, PhiSTLOption::PHI_NEW);
  LBSVecOps::SetPrimarySTLvectorFromMultiGSPETScVec(
    *lbs_problem, groupset_ids, x_, PhiSTLOption::PHI_OLD);

  // Compute final k_eff
  double k_eff = ComputeFissionProduction(*lbs_problem, lbs_problem->GetPhiNewLocal());
  nl_context_ptr->kresid_func_context.k_eff = k_eff;

  PetscInt number_of_func_evals = 0;
  OpenSnPETScCall(SNESGetNumberFunctionEvals(nl_solver_, &number_of_func_evals));
  SNESConvergedReason reason = SNES_CONVERGED_ITERATING;
  OpenSnPETScCall(SNESGetConvergedReason(nl_solver_, &reason));
  const auto status = SNESReasonToPETScSolverStatus(reason);

  // Print summary
  std::stringstream summary;
  summary << FormatKEigenFinalSummary("NLKE",
                                      k_eff,
                                      -1.0,
                                      static_cast<std::size_t>(number_of_func_evals),
                                      "func evals",
                                      PETScSolverStatusToIterationStatus(status));
  if (log.GetVerbosity() >= 2)
    summary << ", detail = " << GetPETScConvergedReasonstring(reason);
  log.Log() << program_timer.GetTimeString() << " " << summary.str();
}

} // namespace opensn
