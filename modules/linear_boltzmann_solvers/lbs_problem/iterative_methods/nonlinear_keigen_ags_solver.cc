// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "modules/linear_boltzmann_solvers/lbs_problem/iterative_methods/nonlinear_keigen_ags_solver.h"
#include "modules/linear_boltzmann_solvers/lbs_problem/iterative_methods/snes_k_monitor.h"
#include "modules/linear_boltzmann_solvers/lbs_problem/iterative_methods/nonlinear_keigen_ags_residual_func.h"
#include "modules/linear_boltzmann_solvers/lbs_problem/lbs_vecops.h"
#include "modules/linear_boltzmann_solvers/lbs_problem/lbs_compute.h"
#include "framework/math/petsc_utils/petsc_utils.h"
#include "framework/runtime.h"
#include "framework/logging/log.h"
#include <petscsnes.h>
#include <iomanip>

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
  for (auto& groupset : lbs_problem->GetGroupsets())
    nl_context_ptr->groupset_ids.push_back(groupset.id);
}

void
NLKEigenvalueAGSSolver::SetMonitor()
{
  auto nl_context_ptr = GetNLKAGSContextPtr(context_ptr_, __PRETTY_FUNCTION__);

  auto& lbs_problem = nl_context_ptr->lbs_problem;
  if (lbs_problem->GetOptions().verbose_outer_iterations)
    SNESMonitorSet(nl_solver_, &KEigenSNESMonitor, &nl_context_ptr->kresid_func_context, nullptr);

  if (lbs_problem->GetOptions().verbose_inner_iterations)
  {
    KSP ksp = nullptr;
    SNESGetKSP(nl_solver_, &ksp);
    KSPMonitorSet(ksp, &KEigenKSPMonitor, &nl_context_ptr->kresid_func_context, nullptr);
  }
}

void
NLKEigenvalueAGSSolver::SetSystemSize()
{
  auto nl_context_ptr = GetNLKAGSContextPtr(context_ptr_, __PRETTY_FUNCTION__);

  auto& lbs_problem = nl_context_ptr->lbs_problem;
  auto sizes = lbs_problem->GetNumPhiIterativeUnknowns();

  num_local_dofs_ = static_cast<int64_t>(sizes.first);
  num_global_dofs_ = static_cast<int64_t>(sizes.second);
}

void
NLKEigenvalueAGSSolver::SetSystem()
{
  // Create the vectors
  x_ = CreateVector(num_local_dofs_, num_global_dofs_);
  VecDuplicate(x_, &r_);
}

void
NLKEigenvalueAGSSolver::SetFunction()
{
  auto nl_context_ptr = GetNLKAGSContextPtr(context_ptr_, __PRETTY_FUNCTION__);

  SNESSetFunction(nl_solver_, r_, NLKEigenResidualFunction, &nl_context_ptr->kresid_func_context);
}

void
NLKEigenvalueAGSSolver::SetJacobian()
{
  MatCreateSNESMF(nl_solver_, &J_);
  SNESSetJacobian(nl_solver_, J_, J_, MatMFFDComputeJacobian, nullptr);
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

  // Unpack solution
  const auto& groups = lbs_problem->GetGroups();
  LBSVecOps::SetPrimarySTLvectorFromGroupScopedPETScVec(
    *lbs_problem, groups.front().id, groups.back().id, x_, PhiSTLOption::PHI_NEW);

  // Compute final k_eff
  double k_eff = ComputeFissionProduction(*lbs_problem, lbs_problem->GetPhiOldLocal());

  PetscInt number_of_func_evals = 0;
  SNESGetNumberFunctionEvals(nl_solver_, &number_of_func_evals);

  // Print summary
  log.Log() << "\n"
            << "        Final k-eigenvalue    :        " << std::fixed << std::setw(10)
            << std::setprecision(7) << k_eff << " (Number of Sweeps:" << number_of_func_evals << ")"
            << "\n";
}

} // namespace opensn
