#include "modules/linear_boltzmann_solvers/a_lbs_solver/acceleration/nl_keigen_acc_solver.h"

#include "modules/linear_boltzmann_solvers/a_lbs_solver/acceleration/nl_keigen_acc_residual_func.h"

#include "modules/linear_boltzmann_solvers/a_lbs_solver/iterative_methods/snes_k_monitor.h"

#include "framework/math/petsc_utils/petsc_utils.h"

#include "framework/runtime.h"
#include "framework/logging/log.h"

#include <iomanip>

namespace opensn
{
namespace lbs
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

  if (nl_context_ptr->verbosity_level_ >= 1)
    SNESMonitorSet(
      nl_solver_, &lbs::KEigenSNESMonitor, &nl_context_ptr->kresid_func_context_, nullptr);

  if (nl_context_ptr->verbosity_level_ >= 2)
  {
    KSP ksp;
    SNESGetKSP(nl_solver_, &ksp);
    KSPMonitorSet(ksp, &lbs::KEigenKSPMonitor, &nl_context_ptr->kresid_func_context_, nullptr);
  }
}

void
NLKEigenDiffSolver::SetSystemSize()
{
  auto nl_context_ptr = GetNLKDiffContextPtr(context_ptr_, __PRETTY_FUNCTION__);

  auto& diff_solver = nl_context_ptr->diff_solver_;
  auto sizes = diff_solver.GetNumPhiIterativeUnknowns();

  num_local_dofs_ = static_cast<int64_t>(sizes.first);
  num_globl_dofs_ = static_cast<int64_t>(sizes.second);
}

void
NLKEigenDiffSolver::SetSystem()
{
  // Create the vectors
  x_ = CreateVector(num_local_dofs_, num_globl_dofs_);
  VecDuplicate(x_, &r_);
}

void
NLKEigenDiffSolver::SetFunction()
{
  auto nl_context_ptr = GetNLKDiffContextPtr(context_ptr_, __PRETTY_FUNCTION__);

  SNESSetFunction(
    nl_solver_, r_, NLKEigenAccResidualFunction, &nl_context_ptr->kresid_func_context_);
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

  auto& lbs_solver = nl_context_ptr->lbs_solver_;
  auto& groupsets = lbs_solver.Groupsets();
  auto& front_gs = groupsets.front();

  auto& phi_old_local = lbs_solver.PhiOldLocal();
  auto& phi_new_local = lbs_solver.PhiNewLocal();

  auto delta_phi = nl_context_ptr->PhiVecToSTLVec(x_);
  auto& phi_lph_ip1 = nl_context_ptr->phi_lph_ip1_;

  auto phi_lp1_temp = phi_lph_ip1 + delta_phi;
  lbs_solver.GSProjectBackPhi0(front_gs, phi_lp1_temp, phi_new_local);
  lbs_solver.GSScopedCopyPrimarySTLvectors(front_gs, phi_new_local, phi_old_local);

  // Compute final k_eff
  double k_eff = nl_context_ptr->kresid_func_context_.k_eff;

  const double production = lbs_solver.ComputeFissionProduction(phi_old_local);
  lbs_solver.ScalePhiVector(PhiSTLOption::PHI_OLD, k_eff / production);

  PetscInt number_of_func_evals;
  SNESGetNumberFunctionEvals(nl_solver_, &number_of_func_evals);

  // Print summary
  if (nl_context_ptr->verbosity_level_ >= 1)
    log.Log() << "        Final lambda-eigenvalue    :        " << std::fixed << std::setw(10)
              << std::setprecision(7) << k_eff << " (num_DOps:" << number_of_func_evals << ")"
              << "\n";
}

} // namespace lbs
} // namespace opensn
