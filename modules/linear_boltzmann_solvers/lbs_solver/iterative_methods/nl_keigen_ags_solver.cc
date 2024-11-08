// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "modules/linear_boltzmann_solvers/lbs_solver/iterative_methods/nl_keigen_ags_solver.h"

#include "modules/linear_boltzmann_solvers/lbs_solver/iterative_methods/snes_k_monitor.h"

#include "modules/linear_boltzmann_solvers/lbs_solver/iterative_methods/nl_keigen_ags_residual_func.h"

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

  auto& lbs_solver = nl_context_ptr->lbs_solver;
  for (auto& groupset : lbs_solver.Groupsets())
    nl_context_ptr->groupset_ids.push_back(groupset.id);
}

void
NLKEigenvalueAGSSolver::SetMonitor()
{
  auto nl_context_ptr = GetNLKAGSContextPtr(context_ptr_, __PRETTY_FUNCTION__);

  auto& lbs_solver = nl_context_ptr->lbs_solver;
  if (lbs_solver.Options().verbose_outer_iterations)
    SNESMonitorSet(nl_solver_, &KEigenSNESMonitor, &nl_context_ptr->kresid_func_context, nullptr);

  if (lbs_solver.Options().verbose_inner_iterations)
  {
    KSP ksp;
    SNESGetKSP(nl_solver_, &ksp);
    KSPMonitorSet(ksp, &KEigenKSPMonitor, &nl_context_ptr->kresid_func_context, nullptr);
  }
}

void
NLKEigenvalueAGSSolver::SetSystemSize()
{
  auto nl_context_ptr = GetNLKAGSContextPtr(context_ptr_, __PRETTY_FUNCTION__);

  auto& lbs_solver = nl_context_ptr->lbs_solver;
  auto sizes = lbs_solver.GetNumPhiIterativeUnknowns();

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

  auto& lbs_solver = nl_context_ptr->lbs_solver;
  const auto& groupset_ids = nl_context_ptr->groupset_ids;

  lbs_solver.SetMultiGSPETScVecFromPrimarySTLvector(groupset_ids, x_, PhiSTLOption::PHI_OLD);
}

void
NLKEigenvalueAGSSolver::PostSolveCallback()
{
  auto nl_context_ptr = GetNLKAGSContextPtr(context_ptr_, __PRETTY_FUNCTION__);

  auto& lbs_solver = nl_context_ptr->lbs_solver;

  // Unpack solution
  const auto& groups = lbs_solver.Groups();
  lbs_solver.SetPrimarySTLvectorFromGroupScopedPETScVec(
    groups.front().id, groups.back().id, x_, lbs_solver.PhiNewLocal());

  // Compute final k_eff
  double k_eff = lbs_solver.ComputeFissionProduction(lbs_solver.PhiNewLocal());

  PetscInt number_of_func_evals;
  SNESGetNumberFunctionEvals(nl_solver_, &number_of_func_evals);

  // Print summary
  log.Log() << "\n"
            << "        Final k-eigenvalue    :        " << std::fixed << std::setw(10)
            << std::setprecision(7) << k_eff << " (Number of Sweeps:" << number_of_func_evals << ")"
            << "\n";
}

} // namespace opensn
