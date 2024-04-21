// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include <petscksp.h>
#include <vector>

namespace opensn
{
class ParallelVector;

/**Generalized solver structure.*/
struct PETScSolverSetup
{
  KSP ksp;
  PC pc;

  std::string solver_name = "KSPSolver";

  std::string solver_type = KSPGMRES;
  std::string preconditioner_type = PCNONE;

  double relative_residual_tol = 1.0e-6;
  int maximum_iterations = 100;
};

/**
 * Creates a general vector.
 *
 * This is a macro for:
 * \code
 * Vec x;
 * VecCreate(opensn::mpi_comm,&x);
 * VecSetType(x,VECMPI);
 * VecSetSizes(x, local_size, global_size);
 * VecSetOption(x,VEC_IGNORE_NEGATIVE_INDICES,PETSC_TRUE);
 *
 * return x;
 * \endcode
 */
Vec CreateVector(int64_t local_size, int64_t global_size);

/**
 * Creates a general vector.
 *
 * This is a function for:
 * \code
 * VecCreate(opensn::mpi_comm,&x);
 * VecSetType(x,VECMPI);
 * VecSetSizes(x, local_size, global_size);
 * VecSetOption(x,VEC_IGNORE_NEGATIVE_INDICES,PETSC_TRUE);
 * \endcode
 */
void CreateVector(Vec& x, int64_t local_size, int64_t global_size);

/**
 * Creates a general vector with ghost value support.
 *
 * This is a macro for:
 * \code
 * Vec x;
 * VecCreateGhost(opensn::mpi_comm,
 *                local_size,
 *                global_size,
 *                nghosts,
 *                ghost_indices.data(),
 *                &x);
 *
 * VecSetOption(x,VEC_IGNORE_NEGATIVE_INDICES,PETSC_TRUE);
 *
 * return x;
 * \endcode
 */
Vec CreateVectorWithGhosts(int64_t local_size,
                           int64_t global_size,
                           int64_t nghosts,
                           const std::vector<int64_t>& ghost_indices);

/**
 * Creates a general square matrix.
 *
 * This is a function for:
 * \code
 * Mat A;
 * MatCreate(opensn::mpi_comm,&A);
 * MatSetType(A,MATMPIAIJ);
 * MatSetSizes(A,local_size, local_size,
 *               global_size, global_size);
 *
 * MatMPIAIJSetPreallocation(A,1, nullptr,
 *                           0, nullptr);
 * MatSetOption(A, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE);
 * MatSetOption(A, MAT_IGNORE_ZERO_ENTRIES, PETSC_TRUE);
 *
 * return A;
 * \endcode
 */
Mat CreateSquareMatrix(int64_t local_size, int64_t global_size);

/**
 * Creates a general square matrix.
 *
 * This is a function for:
 * \code
 * MatCreate(opensn::mpi_comm,&A);
 * MatSetType(A,MATMPIAIJ);
 * MatSetSizes(A,local_size, local_size,
 *               global_size, global_size);
 *
 * MatMPIAIJSetPreallocation(A,1, nullptr,
 *                           0, nullptr);
 * MatSetOption(A, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE);
 * MatSetOption(A, MAT_IGNORE_ZERO_ENTRIES, PETSC_TRUE);
 * \endcode
 */
void CreateSquareMatrix(Mat& A, int64_t local_size, int64_t global_size);

/**
 * Initializes the sparsity pattern of a matrix.
 *
 * This is a macro for:
 * \code
 * MatMPIAIJSetPreallocation(A,0,nodal_nnz_in_diag.data(),
 *                             0,nodal_nnz_off_diag.data());
 * MatSetOption(A, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE);
 * MatSetOption(A, MAT_IGNORE_ZERO_ENTRIES, PETSC_TRUE);
 * MatSetUp(A);
 * \endcode
 */
void InitMatrixSparsity(Mat& A,
                        const std::vector<int64_t>& nodal_nnz_in_diag,
                        const std::vector<int64_t>& nodal_nnz_off_diag);

/**
 * Initializes the sparsity pattern of a matrix.
 *
 * This is a macro for:
 * \code
 * MatMPIAIJSetPreallocation(A,0,nodal_nnz_in_diag.data(),
 *                             0,nodal_nnz_off_diag.data());
 * MatSetOption(A, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE);
 * MatSetOption(A, MAT_IGNORE_ZERO_ENTRIES, PETSC_TRUE);
 * MatSetUp(A);
 * \endcode
 */
void InitMatrixSparsity(Mat& A, int64_t nodal_nnz_in_diag, int64_t nodal_nnz_off_diag);

/**
 * Creates a common Krylov-solver setup.
 *
 * This is a macro for:
 * \code
 * PETScSolverSetup setup;
 *
 * KSPCreate(opensn::mpi_comm,&setup.ksp);
 * KSPSetOperators(setup.ksp,ref_matrix,ref_matrix);
 * KSPSetType(setup.ksp,in_solver_type.c_str());
 *
 * KSPSetOptionsPrefix(setup.ksp,in_solver_name.c_str());
 *
 * KSPGetPC(setup.ksp,&setup.pc);
 * PCSetType(setup.pc,in_preconditioner_type.c_str());
 *
 * KSPSetTolerances(setup.ksp,1.e-50,
 *                  in_relative_residual_tolerance,1.0e50,
 *                  in_maximum_iterations);
 * KSPSetInitialGuessNonzero(setup.ksp,PETSC_TRUE);
 *
 * KSPMonitorSet(setup.ksp,&KSPMonitorRelativeToRHS,NULL,NULL);
 *
 * return setup;
 * \endcode
 */
PETScSolverSetup CreateCommonKrylovSolverSetup(Mat matrix,
                                               const std::string& solver_name = "KSPSolver",
                                               const std::string& solver_type = KSPGMRES,
                                               const std::string& preconditioner_type = PCNONE,
                                               double rel_tol = PETSC_DEFAULT,
                                               double abs_tol = PETSC_DEFAULT,
                                               int64_t maximum_iterations = 100);

/**
 * General monitor that print the residual norm relative to the right-hand side norm.
 */
PetscErrorCode KSPMonitorRelativeToRHS(KSP ksp, PetscInt n, PetscReal rnorm, void*);

/**
 * Copies a PETSc vector to a STL vector. Only the local portion is copied.
 */
void CopyVecToSTLvector(Vec x, std::vector<double>& data, size_t N, bool resize_STL = true);

/**
 * Copies a PETSc vector to a STL vector. Only the local portion is copied.
 */
void
CopyVecToSTLvectorWithGhosts(Vec x, std::vector<double>& data, size_t N, bool resize_STL = true);

void CopySTLvectorToVec(const std::vector<double>& data, Vec x, size_t N);
void CopyParallelVectorToVec(const ParallelVector& y, Vec x);

/**
 * Copies global values from a PETSc vector to a STL vector.
 */
void CopyGlobalVecToSTLvector(Vec x,
                              const std::vector<int64_t>& global_indices,
                              std::vector<double>& data);

/**
 * Communicates ghost entries of a ghost vector. This operation is suitable when only a single
 * vector is communicated. When more than vector is communicated it would be more efficient to
 * "Begin" all the vectors followed by and "End" of each vector.
 */
void CommunicateGhostEntries(Vec x);

/**Simple data structure to keep track of a ghost vector's
 * localized views.*/
struct GhostVecLocalRaw
{
  Vec x_localized = nullptr;
  double* x_localized_raw = nullptr;

  /**Returns a copy of the value at the specified index.*/
  double operator[](int index) { return x_localized_raw[index]; }

  /**Returns a reference of the value at the specified index.*/
  double& operator()(int index) { return x_localized_raw[index]; }
};

/**
 * Gets a local raw view of a ghost vector.
 */
GhostVecLocalRaw GetGhostVectorLocalViewRead(Vec x);

/**
 * Gets a local raw view of a ghost vector.
 */
void RestoreGhostVectorLocalViewRead(Vec x, GhostVecLocalRaw& local_data);

/**
 * Gets the string value of a converged reason.
 */
std::string GetPETScConvergedReasonstring(KSPConvergedReason reason);

} // namespace opensn
