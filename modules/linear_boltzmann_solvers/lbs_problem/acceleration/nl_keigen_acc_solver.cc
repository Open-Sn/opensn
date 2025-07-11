// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "modules/linear_boltzmann_solvers/lbs_problem/acceleration/nl_keigen_acc_solver.h"
#include "modules/linear_boltzmann_solvers/lbs_problem/acceleration/nl_keigen_acc_residual_func.h"
#include "modules/linear_boltzmann_solvers/lbs_problem/iterative_methods/snes_k_monitor.h"
#include "modules/linear_boltzmann_solvers/lbs_problem/lbs_vecops.h"
#include "modules/linear_boltzmann_solvers/lbs_problem/lbs_compute.h"
#include "modules/linear_boltzmann_solvers/lbs_problem/acceleration/wgdsa.h"
#include "framework/math/petsc_utils/petsc_utils.h"
#include "framework/runtime.h"
#include "framework/logging/log.h"
#include <iomanip>

namespace opensn
{

namespace
{

std::shared_ptr<NLKEigenDiffContext>
GetNLKDiffContextPtr(const std::shared_ptr<NonLinearSolverContext>& context,
                     const std::string& func_name)
{
  auto nlk_eigen_diff_context = std::dynamic_pointer_cast<NLKEigenDiffContext>(context);
  if (not nlk_eigen_diff_context)
    throw std::runtime_error(std::string(func_name) + ": context casting failure");
  return nlk_eigen_diff_context;
}

} // namespace

void
NLKEigenDiffSolver::SetMonitor()
{
  auto nl_context_ptr = GetNLKDiffContextPtr(context_ptr_, __PRETTY_FUNCTION__);

  if (nl_context_ptr->verbosity_level >= 1)
    SNESMonitorSet(nl_solver_, &KEigenSNESMonitor, &nl_context_ptr->kresid_func_context, nullptr);

  if (nl_context_ptr->verbosity_level >= 2)
  {
    KSP ksp;
    SNESGetKSP(nl_solver_, &ksp);
    KSPMonitorSet(ksp, &KEigenKSPMonitor, &nl_context_ptr->kresid_func_context, nullptr);
  }
}

void
NLKEigenDiffSolver::SetSystemSize()
{
  auto nl_context_ptr = GetNLKDiffContextPtr(context_ptr_, __PRETTY_FUNCTION__);

  auto& diff_solver = nl_context_ptr->diff_solver;
  auto sizes = diff_solver.GetNumPhiIterativeUnknowns();

  num_local_dofs_ = static_cast<int64_t>(sizes.first);
  num_global_dofs_ = static_cast<int64_t>(sizes.second);
}

void
NLKEigenDiffSolver::SetSystem()
{
  // Create the vectors
  x_ = CreateVector(num_local_dofs_, num_global_dofs_);
  VecDuplicate(x_, &r_);
}

void
NLKEigenDiffSolver::SetFunction()
{
  auto nl_context_ptr = GetNLKDiffContextPtr(context_ptr_, __PRETTY_FUNCTION__);

  SNESSetFunction(
    nl_solver_, r_, NLKEigenAccResidualFunction, &nl_context_ptr->kresid_func_context);
}

void
NLKEigenDiffSolver::SetJacobian()
{
  MatCreateSNESMF(nl_solver_, &J_);
  SNESSetJacobian(nl_solver_, J_, J_, MatMFFDComputeJacobian, nullptr);
}

void
NLKEigenDiffSolver::SetInitialGuess()
{
  VecSet(x_, 0.0);
}

void
NLKEigenDiffSolver::PostSolveCallback()
{
  auto nl_context_ptr = GetNLKDiffContextPtr(context_ptr_, __PRETTY_FUNCTION__);

  auto& lbs_problem = nl_context_ptr->lbs_problem;
  auto& groupsets = lbs_problem.GetGroupsets();
  auto& front_gs = groupsets.front();

  auto& phi_old_local = lbs_problem.GetPhiOldLocal();
  auto& phi_new_local = lbs_problem.GetPhiNewLocal();

  auto delta_phi = nl_context_ptr->PhiVecToSTLVec(x_);
  auto& phi_lph_ip1 = nl_context_ptr->phi_lph_ip1;

  auto phi_lp1_temp = phi_lph_ip1 + delta_phi;
  WGDSA::GSProjectBackPhi0(lbs_problem, front_gs, phi_lp1_temp, phi_new_local);
  LBSVecOps::GSScopedCopyPrimarySTLvectors(lbs_problem, front_gs, phi_new_local, phi_old_local);

  // Compute final k_eff
  double k_eff = nl_context_ptr->kresid_func_context.k_eff;

  const double production = ComputeFissionProduction(lbs_problem, phi_old_local);
  LBSVecOps::ScalePhiVector(lbs_problem, PhiSTLOption::PHI_OLD, k_eff / production);

  PetscInt number_of_func_evals;
  SNESGetNumberFunctionEvals(nl_solver_, &number_of_func_evals);

  // Print summary
  if (nl_context_ptr->verbosity_level >= 1)
    log.Log() << "        Final lambda-eigenvalue    :        " << std::fixed << std::setw(10)
              << std::setprecision(7) << k_eff << " (num_DOps:" << number_of_func_evals << ")"
              << "\n";
}

} // namespace opensn
