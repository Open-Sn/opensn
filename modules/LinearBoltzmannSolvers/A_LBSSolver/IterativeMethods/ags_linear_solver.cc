#include "ags_linear_solver.h"

#include "opensn/modules/LinearBoltzmannSolvers/A_LBSSolver/lbs_solver.h"

#include "opensn/framework/math/PETScUtils/petsc_utils.h"
#include "opensn/framework/math/LinearSolver/linear_matrix_action_Ax.h"

#include <petscksp.h>

#include "opensn/framework/chi_runtime.h"
#include "opensn/framework/logging/chi_log.h"

#include <iomanip>

#define sc_int64_t static_cast<int64_t>

#define GetAGSContextPtr(x) std::dynamic_pointer_cast<AGSContext<Mat, Vec, KSP>>(x)
namespace lbs
{

template <>
void
AGSLinearSolver<Mat, Vec, KSP>::SetSystemSize()
{
  auto ags_context_ptr = GetAGSContextPtr(context_ptr_);

  const auto sizes = ags_context_ptr->SystemSize();

  num_local_dofs_ = sizes.first;
  num_globl_dofs_ = sizes.second;
}

template <>
void
AGSLinearSolver<Mat, Vec, KSP>::SetSystem()
{
  x_ = chi_math::PETScUtils::CreateVector(sc_int64_t(num_local_dofs_), sc_int64_t(num_globl_dofs_));

  VecSet(x_, 0.0);
  VecDuplicate(x_, &b_);

  //============================================= Create the matrix-shell
  MatCreateShell(PETSC_COMM_WORLD,
                 sc_int64_t(num_local_dofs_),
                 sc_int64_t(num_local_dofs_),
                 sc_int64_t(num_globl_dofs_),
                 sc_int64_t(num_globl_dofs_),
                 &(*context_ptr_),
                 &A_);

  //============================================= Set the action-operator
  MatShellSetOperation(A_, MATOP_MULT, (void (*)())chi_math::LinearSolverMatrixAction<Mat, Vec>);

  //============================================= Set solver operators
  KSPSetOperators(solver_, A_, A_);
  KSPSetUp(solver_);
}

template <>
void
AGSLinearSolver<Mat, Vec, KSP>::SetPreconditioner()
{
  auto ags_context_ptr = GetAGSContextPtr(context_ptr_);

  ags_context_ptr->SetPreconditioner(solver_);
}

template <>
void
AGSLinearSolver<Mat, Vec, KSP>::SetRHS()
{
}

template <>
void
AGSLinearSolver<Mat, Vec, KSP>::SetInitialGuess()
{
}

template <>
void
AGSLinearSolver<Mat, Vec, KSP>::Solve()
{
  auto ags_context_ptr = GetAGSContextPtr(context_ptr_);
  auto& lbs_solver = ags_context_ptr->lbs_solver_;

  const int gid_i = GroupSpanFirstID();
  const int gid_f = GroupSpanLastID();
  const auto& phi = lbs_solver.PhiOldLocal();

  Vec x_old;
  VecDuplicate(x_, &x_old);

  // Save qmoms to be restored after each iteration.
  // This is necessary for multiple ags iterations to function
  // and for keigen-value problems
  const auto saved_qmoms = lbs_solver.QMomentsLocal();

  for (int iter = 0; iter < tolerance_options_.maximum_iterations; ++iter)
  {

    lbs_solver.SetGroupScopedPETScVecFromPrimarySTLvector(gid_i, gid_f, x_old, phi);

    for (auto& solver : ags_context_ptr->sub_solvers_list_)
    {
      solver->Setup();
      solver->Solve();
    }

    lbs_solver.SetGroupScopedPETScVecFromPrimarySTLvector(gid_i, gid_f, x_, phi);

    VecAXPY(x_old, -1.0, x_);
    PetscReal error_norm;
    VecNorm(x_old, NORM_2, &error_norm);
    PetscReal sol_norm;
    VecNorm(x_, NORM_2, &sol_norm);

    if (verbose_)
      Chi::log.Log() << "********** AGS solver iteration " << std::setw(3) << iter << " "
                     << " Relative change " << std::setw(10) << std::setprecision(4)
                     << error_norm / sol_norm;

    lbs_solver.QMomentsLocal() = saved_qmoms; // Restore qmoms

    if (error_norm < tolerance_options_.residual_absolute) break;
  } // for iteration

  VecDestroy(&x_old);
}

template <>
AGSLinearSolver<Mat, Vec, KSP>::~AGSLinearSolver()
{
  MatDestroy(&A_);
}

} // namespace lbs
