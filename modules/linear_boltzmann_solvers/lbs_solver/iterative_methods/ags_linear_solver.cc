// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "modules//linear_boltzmann_solvers/lbs_solver/iterative_methods/ags_linear_solver.h"
#include "modules/linear_boltzmann_solvers/lbs_solver/lbs_solver.h"
#include "framework/math/linear_solver/linear_matrix_action_Ax.h"
#include "framework/math/petsc_utils/petsc_utils.h"
#include "framework/logging/log.h"
#include "framework/runtime.h"
#include "caliper/cali.h"
#include <petscksp.h>
#include <iomanip>

namespace opensn
{
namespace lbs
{

void
AGSLinearSolver::SetSystemSize()
{
  CALI_CXX_MARK_SCOPE("AGSLinearSolver::SetSystemSize");

  auto ags_context_ptr = std::dynamic_pointer_cast<AGSContext>(context_ptr_);

  const auto sizes = ags_context_ptr->SystemSize();

  num_local_dofs_ = sizes.first;
  num_global_dofs_ = sizes.second;
}

void
AGSLinearSolver::SetSystem()
{
  CALI_CXX_MARK_SCOPE("AGSLinearSolver::SetSystem");

  x_ = CreateVector(static_cast<int64_t>(num_local_dofs_), static_cast<int64_t>(num_global_dofs_));

  VecSet(x_, 0.0);
  VecDuplicate(x_, &b_);

  // Create the matrix-shell
  MatCreateShell(opensn::mpi_comm,
                 static_cast<int64_t>(num_local_dofs_),
                 static_cast<int64_t>(num_local_dofs_),
                 static_cast<int64_t>(num_global_dofs_),
                 static_cast<int64_t>(num_global_dofs_),
                 &(*context_ptr_),
                 &A_);

  // Set the action-operator
  MatShellSetOperation(A_, MATOP_MULT, (void (*)())LinearSolverMatrixAction);

  // Set solver operators
  KSPSetOperators(ksp_, A_, A_);
  KSPSetUp(ksp_);
}

void
AGSLinearSolver::SetPreconditioner()
{
  CALI_CXX_MARK_SCOPE("AGSLinearSolver::SetPreconditioner");

  auto ags_context_ptr = std::dynamic_pointer_cast<AGSContext>(context_ptr_);
  ags_context_ptr->SetPreconditioner(ksp_);
}

void
AGSLinearSolver::SetRHS()
{
}

void
AGSLinearSolver::SetInitialGuess()
{
}

void
AGSLinearSolver::Solve()
{
  CALI_CXX_MARK_SCOPE("AGSLinearSolver::Solve");

  auto ags_context_ptr = std::dynamic_pointer_cast<AGSContext>(context_ptr_);
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
      log.Log() << "********** AGS solver iteration " << std::setw(3) << iter << " "
                << " Relative change " << std::setw(10) << std::setprecision(4)
                << error_norm / sol_norm;

    lbs_solver.QMomentsLocal() = saved_qmoms; // Restore qmoms

    // Write restart data
    if (lbs_solver.RestartsEnabled() and lbs_solver.TriggerRestartDump() and
        lbs_solver.Options().enable_ags_restart_write)
    {
      lbs_solver.WriteRestartData();
    }

    if (error_norm < tolerance_options_.residual_absolute)
      break;
  } // for iteration

  // If restarts are enabled, always write a restart dump upon convergence or
  // when we reach the iteration limit
  if (lbs_solver.RestartsEnabled() && lbs_solver.Options().enable_ags_restart_write)
    lbs_solver.WriteRestartData();

  VecDestroy(&x_old);
}

AGSLinearSolver::~AGSLinearSolver()
{
  MatDestroy(&A_);
}

} // namespace lbs
} // namespace opensn
