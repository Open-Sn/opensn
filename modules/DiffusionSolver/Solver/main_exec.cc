#include "diffusion_solver.h"

#include "mesh/MeshContinuum/chi_meshcontinuum.h"

#include "chi_runtime.h"
#include "chi_log.h"

#include "utils/chi_timer.h"

//###################################################################
/**Executes the diffusion solver using the PETSc library.*/
int
chi_diffusion::Solver::ExecuteS(bool suppress_assembly, bool suppress_solve)
{
  t_assembly_.Reset();

  if (Chi::material_stack.empty())
  {
    Chi::log.Log0Error() << "No materials added to simulation. Add materials.";
    exit(0);
  }

  VecSet(x_, 0.0);
  VecSet(b_, 0.0);

  if (!suppress_assembly)
    Chi::log.Log() << Chi::program_timer.GetTimeString() << " " << TextName()
                   << ": Assembling A locally";

  //================================================== Loop over locally owned
  //                                                   cells
  auto fem_method = basic_options_("discretization_method").StringValue();
  if (fem_method == "PWLC")
  {
    if (!suppress_assembly)
      for (auto& cell : grid_ptr_->local_cells)
        CFEM_Assemble_A_and_b(cell, gi_);
    else {}
  }
  else if (fem_method == "PWLD_MIP")
  {
    if (!suppress_assembly)
      for (auto& cell : grid_ptr_->local_cells)
        PWLD_Assemble_A_and_b(cell, gi_);
    else
      for (auto& cell : grid_ptr_->local_cells)
        PWLD_Assemble_b(cell, gi_);
  }
  else
  {
    Chi::log.Log() << "Diffusion Solver: Finite Element Discretization "
                      "method not specified.";
    Chi::Exit(EXIT_FAILURE);
  }

  if (!suppress_assembly)
    Chi::log.Log() << Chi::program_timer.GetTimeString() << " " << TextName()
                   << ": Done Assembling A locally";
  MPI_Barrier(Chi::mpi.comm);

  //=================================== Call matrix assembly
  if (verbose_info_ || Chi::log.GetVerbosity() >= chi::ChiLog::LOG_0VERBOSE_1)
    Chi::log.Log() << Chi::program_timer.GetTimeString() << " " << TextName()
                   << ": Communicating matrix assembly";

  if (!suppress_assembly)
  {
    Chi::log.Log() << Chi::program_timer.GetTimeString() << " " << TextName()
                   << ": Assembling A globally";
    MatAssemblyBegin(A_, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(A_, MAT_FINAL_ASSEMBLY);

    //================================= Matrix symmetry check
    //    PetscBool is_symmetric;
    //    ierr = MatIsSymmetric(A,1.0e-4,&is_symmetric);
    //    if (!is_symmetric)
    //    {
    //      chi::log.Log0Warning()
    //        << "Assembled matrix is not symmetric";
    //    }

    //================================= Matrix diagonal check
    Chi::log.Log() << Chi::program_timer.GetTimeString() << " " << TextName() << ": Diagonal check";
    PetscBool missing_diagonal;
    PetscInt row;
    MatMissingDiagonal(A_, &missing_diagonal, &row);
    if (missing_diagonal)
      Chi::log.LogAllError() << Chi::program_timer.GetTimeString() << " " << TextName()
                             << ": Missing diagonal detected";

    //================================= Matrix sparsity info
    MatInfo info;
    ierr_ = MatGetInfo(A_, MAT_GLOBAL_SUM, &info);

    Chi::log.Log() << "Number of mallocs used = " << info.mallocs
                   << "\nNumber of non-zeros allocated = " << info.nz_allocated
                   << "\nNumber of non-zeros used = " << info.nz_used
                   << "\nNumber of unneeded non-zeros = " << info.nz_unneeded;
  }
  if (verbose_info_ || Chi::log.GetVerbosity() >= chi::ChiLog::LOG_0VERBOSE_1)
    Chi::log.Log() << Chi::program_timer.GetTimeString() << " " << TextName()
                   << ": Assembling x and b";
  VecAssemblyBegin(x_);
  VecAssemblyEnd(x_);
  VecAssemblyBegin(b_);
  VecAssemblyEnd(b_);

  time_assembly_ = t_assembly_.GetTime() / 1000.0;

  //=================================== Execute solve
  if (suppress_solve)
  {
    Chi::log.Log() << Chi::program_timer.GetTimeString() << " " << TextName()
                   << ": Setting up solver and preconditioner\n";
    PCSetUp(pc_);
    KSPSetUp(ksp_);
  }
  else
  {
    if (verbose_info_ || Chi::log.GetVerbosity() >= chi::ChiLog::LOG_0VERBOSE_1)
      Chi::log.Log() << Chi::program_timer.GetTimeString() << " " << TextName()
                     << ": Solving system\n";
    t_solve_.Reset();
    PCSetUp(pc_);
    KSPSetUp(ksp_);
    KSPSolve(ksp_, b_, x_);
    time_solve_ = t_solve_.GetTime() / 1000.0;

    //=================================== Populate field vector
    //    if (fem_method == "PWLD_MIP" or fem_method == "PWLD_MIP_GAGG")
    //    {
    //      const double* x_ref;
    //      VecGetArrayRead(x,&x_ref);
    //
    //      for (int i=0; i < local_dof_count; i++)
    //        pwld_phi_local[i] = x_ref[i];
    //
    //      VecRestoreArrayRead(x,&x_ref);
    //    }

    //=================================== Get convergence reason
    KSPConvergedReason reason;
    KSPGetConvergedReason(ksp_, &reason);
    if (verbose_info_ || reason != KSP_CONVERGED_RTOL)
      Chi::log.Log() << "Convergence reason: "
                     << chi_physics::GetPETScConvergedReasonstring(reason);

    //=================================== Location wise view
    if (Chi::mpi.location_id == 0)
    {
      int64_t its;
      ierr_ = KSPGetIterationNumber(ksp_, &its);
      Chi::log.Log() << Chi::program_timer.GetTimeString() << " " << TextName() << "[g=" << gi_
                     << "-" << gi_ + G_ - 1 << "]: Number of iterations =" << its;

      if (verbose_info_ || Chi::log.GetVerbosity() >= chi::ChiLog::LOG_0VERBOSE_1)
      {
        Chi::log.Log() << "Timing:";
        Chi::log.Log() << "Assembling the matrix: " << time_assembly_;
        Chi::log.Log() << "Solving the system   : " << time_solve_;
      }
    }

    UpdateFieldFunctions();

    if (verbose_info_ || Chi::log.GetVerbosity() >= chi::ChiLog::LOG_0VERBOSE_1)
      Chi::log.Log() << "Diffusion Solver execution completed!\n";
  } // if not suppressed solve

  return 0;
}