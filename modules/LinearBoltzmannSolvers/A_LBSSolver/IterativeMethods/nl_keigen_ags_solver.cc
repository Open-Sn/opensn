#include "modules/LinearBoltzmannSolvers/A_LBSSolver/IterativeMethods/nl_keigen_ags_solver.h"

#include "modules/LinearBoltzmannSolvers/A_LBSSolver/IterativeMethods/snes_k_monitor.h"

#include "modules/LinearBoltzmannSolvers/A_LBSSolver/IterativeMethods/nl_keigen_ags_residual_func.h"

#include "framework/math/petsc_utils/petsc_utils.h"

#include "framework/runtime.h"
#include "framework/logging/log.h"

#include <petscsnes.h>

#include <iomanip>

#define SNESTypes Mat, Vec, SNES
#define CheckContext(x)                                                                            \
  if (not x)                                                                                       \
  throw std::runtime_error(std::string(__PRETTY_FUNCTION__) + ": context casting failure")
#define GetNLKAGSContextPtr(x)                                                                     \
  std::dynamic_pointer_cast<NLKEigenAGSContext<Vec, SNES>>(x);                                     \
  CheckContext(x)

namespace lbs
{

template <>
void
NLKEigenvalueAGSSolver<SNESTypes>::PreSetupCallback()
{
  auto nl_context_ptr = GetNLKAGSContextPtr(context_ptr_);

  auto& lbs_solver = nl_context_ptr->lbs_solver_;
  for (auto& groupset : lbs_solver.Groupsets())
    nl_context_ptr->groupset_ids.push_back(groupset.id_);
}

template <>
void
NLKEigenvalueAGSSolver<SNESTypes>::SetMonitor()
{
  auto nl_context_ptr = GetNLKAGSContextPtr(context_ptr_);

  auto& lbs_solver = nl_context_ptr->lbs_solver_;
  if (lbs_solver.Options().verbose_outer_iterations)
    SNESMonitorSet(
      nl_solver_, &lbs::KEigenSNESMonitor, &nl_context_ptr->kresid_func_context_, nullptr);

  if (lbs_solver.Options().verbose_inner_iterations)
  {
    KSP ksp;
    SNESGetKSP(nl_solver_, &ksp);
    KSPMonitorSet(ksp, &lbs::KEigenKSPMonitor, &nl_context_ptr->kresid_func_context_, nullptr);
  }
}

template <>
void
NLKEigenvalueAGSSolver<SNESTypes>::SetSystemSize()
{
  auto nl_context_ptr = GetNLKAGSContextPtr(context_ptr_);

  auto& lbs_solver = nl_context_ptr->lbs_solver_;
  auto sizes = lbs_solver.GetNumPhiIterativeUnknowns();

  num_local_dofs_ = static_cast<int64_t>(sizes.first);
  num_globl_dofs_ = static_cast<int64_t>(sizes.second);
}

template <>
void
NLKEigenvalueAGSSolver<SNESTypes>::SetSystem()
{
  // Create the vectors
  x_ = chi_math::PETScUtils::CreateVector(num_local_dofs_, num_globl_dofs_);
  VecDuplicate(x_, &r_);
}

template <>
void
NLKEigenvalueAGSSolver<SNESTypes>::SetFunction()
{
  auto nl_context_ptr = GetNLKAGSContextPtr(context_ptr_);

  SNESSetFunction(nl_solver_, r_, NLKEigenResidualFunction, &nl_context_ptr->kresid_func_context_);
}

template <>
void
NLKEigenvalueAGSSolver<SNESTypes>::SetJacobian()
{
  MatCreateSNESMF(nl_solver_, &J_);
  SNESSetJacobian(nl_solver_, J_, J_, MatMFFDComputeJacobian, nullptr);
}

template <>
void
NLKEigenvalueAGSSolver<SNESTypes>::SetInitialGuess()
{
  auto nl_context_ptr = GetNLKAGSContextPtr(context_ptr_);

  auto& lbs_solver = nl_context_ptr->lbs_solver_;
  const auto& groupset_ids = nl_context_ptr->groupset_ids;

  lbs_solver.SetMultiGSPETScVecFromPrimarySTLvector(groupset_ids, x_, PhiSTLOption::PHI_OLD);
}

template <>
void
NLKEigenvalueAGSSolver<SNESTypes>::PostSolveCallback()
{
  auto nl_context_ptr = GetNLKAGSContextPtr(context_ptr_);

  auto& lbs_solver = nl_context_ptr->lbs_solver_;

  // Unpack solution
  const auto& groups = lbs_solver.Groups();
  lbs_solver.SetPrimarySTLvectorFromGroupScopedPETScVec(
    groups.front().id_, groups.back().id_, x_, lbs_solver.PhiOldLocal());

  // Compute final k_eff
  double k_eff = lbs_solver.ComputeFissionProduction(lbs_solver.PhiOldLocal());

  PetscInt number_of_func_evals;
  SNESGetNumberFunctionEvals(nl_solver_, &number_of_func_evals);

  // Print summary
  Chi::log.Log() << "\n"
                 << "        Final k-eigenvalue    :        " << std::fixed << std::setw(10)
                 << std::setprecision(7) << k_eff << " (num_TrOps:" << number_of_func_evals << ")"
                 << "\n";
}

} // namespace lbs
