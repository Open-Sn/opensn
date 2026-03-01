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
    OpenSnPETScCall(SNESMonitorSet(
      nl_solver_, &KEigenSNESMonitor, &nl_context_ptr->kresid_func_context, nullptr));

  if (nl_context_ptr->verbosity_level >= 2)
  {
    KSP ksp = nullptr;
    OpenSnPETScCall(SNESGetKSP(nl_solver_, &ksp));
    OpenSnPETScCall(
      KSPMonitorSet(ksp, &KEigenKSPMonitor, &nl_context_ptr->kresid_func_context, nullptr));
  }
}

void
NLKEigenDiffSolver::SetSystemSize()
{
  auto nl_context_ptr = GetNLKDiffContextPtr(context_ptr_, __PRETTY_FUNCTION__);

  auto& diff_solver = nl_context_ptr->diff_solver;
  std::tie(num_local_dofs_, num_global_dofs_) = diff_solver.GetNumPhiIterativeUnknowns();
}

void
NLKEigenDiffSolver::SetSystem()
{
  // Create the vectors
  x_ = CreateVector(num_local_dofs_, num_global_dofs_);
  OpenSnPETScCall(VecDuplicate(x_, &r_));
}

void
NLKEigenDiffSolver::SetFunction()
{
  auto nl_context_ptr = GetNLKDiffContextPtr(context_ptr_, __PRETTY_FUNCTION__);

  OpenSnPETScCall(SNESSetFunction(
    nl_solver_, r_, NLKEigenAccResidualFunction, &nl_context_ptr->kresid_func_context));
}

void
NLKEigenDiffSolver::SetJacobian()
{
  OpenSnPETScCall(MatCreateSNESMF(nl_solver_, &J_));
  OpenSnPETScCall(SNESSetJacobian(nl_solver_, J_, J_, MatMFFDComputeJacobian, nullptr));
}

void
NLKEigenDiffSolver::SetInitialGuess()
{
  OpenSnPETScCall(VecSet(x_, 0.0));
}

void
NLKEigenDiffSolver::PostSolveCallback()
{
  auto nl_context_ptr = GetNLKDiffContextPtr(context_ptr_, __PRETTY_FUNCTION__);

  auto& do_problem = nl_context_ptr->do_problem;
  auto& groupsets = do_problem.GetGroupsets();
  auto& front_gs = groupsets.front();

  auto& phi_old_local = do_problem.GetPhiOldLocal();
  auto& phi_new_local = do_problem.GetPhiNewLocal();

  auto delta_phi = nl_context_ptr->PhiVecToSTLVec(x_);
  auto& phi_lph_ip1 = nl_context_ptr->phi_lph_ip1;

  auto phi_lp1_temp = phi_lph_ip1 + delta_phi;
  WGDSA::GSProjectBackPhi0(do_problem, front_gs, phi_lp1_temp, phi_new_local);
  LBSVecOps::GSScopedCopyPrimarySTLvectors(do_problem, front_gs, phi_new_local, phi_old_local);

  // Compute final k_eff
  double k_eff = nl_context_ptr->kresid_func_context.k_eff;

  const double production = ComputeFissionProduction(do_problem, phi_old_local);
  LBSVecOps::ScalePhiVector(do_problem, PhiSTLOption::PHI_OLD, k_eff / production);

  PetscInt number_of_func_evals = 0;
  OpenSnPETScCall(SNESGetNumberFunctionEvals(nl_solver_, &number_of_func_evals));

  // Print summary
  if (nl_context_ptr->verbosity_level >= 1)
    log.Log() << "        Final lambda-eigenvalue    :        " << std::fixed << std::setw(10)
              << std::setprecision(7) << k_eff << " (num_DOps:" << number_of_func_evals << ")"
              << "\n";
}

} // namespace opensn
