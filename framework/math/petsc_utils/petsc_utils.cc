// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "framework/math/petsc_utils/petsc_utils.h"
#include "framework/data_types/parallel_vector/parallel_vector.h"
#include "framework/runtime.h"
#include "framework/logging/log.h"
#include <iomanip>

namespace opensn
{

Vec
CreateVector(int64_t local_size, int64_t global_size)
{
  Vec x = nullptr;
  VecCreate(opensn::mpi_comm, &x);
  VecSetType(x, VECMPI);
  VecSetSizes(x, local_size, global_size);
  VecSetOption(x, VEC_IGNORE_NEGATIVE_INDICES, PETSC_TRUE);

  return x;
}

void
CreateVector(Vec& x, int64_t local_size, int64_t global_size)
{
  VecCreate(opensn::mpi_comm, &x);
  VecSetType(x, VECMPI);
  VecSetSizes(x, local_size, global_size);
  VecSetOption(x, VEC_IGNORE_NEGATIVE_INDICES, PETSC_TRUE);
}

Vec
CreateVectorWithGhosts(int64_t local_size,
                       int64_t global_size,
                       int64_t nghosts,
                       const std::vector<PetscInt>& ghost_indices)
{
  Vec x = nullptr;
  VecCreateGhost(opensn::mpi_comm,
                 local_size,
                 global_size,
                 nghosts,
                 (ghost_indices.empty()) ? NULL : ghost_indices.data(),
                 &x);

  VecSetOption(x, VEC_IGNORE_NEGATIVE_INDICES, PETSC_TRUE);

  return x;
}

Mat
CreateSquareMatrix(int64_t local_size, int64_t global_size)
{
  Mat A = nullptr;
  MatCreate(opensn::mpi_comm, &A);
  MatSetType(A, MATMPIAIJ);
  MatSetSizes(A, local_size, local_size, global_size, global_size);

  MatMPIAIJSetPreallocation(A, 1, nullptr, 0, nullptr);
  MatSetOption(A, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE);
  MatSetOption(A, MAT_IGNORE_ZERO_ENTRIES, PETSC_TRUE);

  return A;
}

void
CreateSquareMatrix(Mat& A, int64_t local_size, int64_t global_size)
{
  MatCreate(opensn::mpi_comm, &A);
  MatSetType(A, MATMPIAIJ);
  MatSetSizes(A, local_size, local_size, global_size, global_size);

  MatMPIAIJSetPreallocation(A, 1, nullptr, 0, nullptr);
  MatSetOption(A, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE);
  MatSetOption(A, MAT_IGNORE_ZERO_ENTRIES, PETSC_TRUE);
}

void
InitMatrixSparsity(Mat& A,
                   const std::vector<int64_t>& nodal_nnz_in_diag,
                   const std::vector<int64_t>& nodal_nnz_off_diag)
{
  MatMPIAIJSetPreallocation(A, 0, nodal_nnz_in_diag.data(), 0, nodal_nnz_off_diag.data());
  MatSetOption(A, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE);
  MatSetOption(A, MAT_IGNORE_ZERO_ENTRIES, PETSC_TRUE);
  MatSetUp(A);
}

void
InitMatrixSparsity(Mat& A, int64_t nodal_nnz_in_diag, int64_t nodal_nnz_off_diag)
{
  MatMPIAIJSetPreallocation(A, nodal_nnz_in_diag, nullptr, nodal_nnz_off_diag, nullptr);
  MatSetOption(A, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE);
  MatSetOption(A, MAT_IGNORE_ZERO_ENTRIES, PETSC_TRUE);
}

PETScSolverSetup
CreateCommonKrylovSolverSetup(Mat matrix,
                              const std::string& solver_name,
                              const std::string& solver_type,
                              const std::string& preconditioner_type,
                              double rel_tol,
                              double abs_tol,
                              int64_t maximum_iterations)
{
  PETScSolverSetup setup;

  KSPCreate(opensn::mpi_comm, &setup.ksp);
  KSPSetOperators(setup.ksp, matrix, matrix);
  KSPSetType(setup.ksp, solver_type.c_str());

  KSPSetOptionsPrefix(setup.ksp, solver_name.c_str());

  KSPGetPC(setup.ksp, &setup.pc);
  PCSetType(setup.pc, preconditioner_type.c_str());

  KSPSetTolerances(setup.ksp, rel_tol, abs_tol, 1.0e50, maximum_iterations);
  KSPSetInitialGuessNonzero(setup.ksp, PETSC_TRUE);

  KSPSetFromOptions(setup.ksp);

  KSPMonitorSet(setup.ksp, &KSPMonitorRelativeToRHS, nullptr, nullptr);

  return setup;
}

PetscErrorCode
KSPMonitorRelativeToRHS(KSP ksp, PetscInt n, PetscReal rnorm, void* /* context */)
{
  Vec Rhs = nullptr;
  KSPGetRhs(ksp, &Rhs);
  double rhs_norm = 0.0;
  VecNorm(Rhs, NORM_2, &rhs_norm);
  if (rhs_norm < 1.0e-12)
    rhs_norm = 1.0;

  // Get solver name
  const char* ksp_name = nullptr;
  KSPGetOptionsPrefix(ksp, &ksp_name);

  // Default to this if ksp_name is NULL
  const char NONAME_SOLVER[] = "NoName-Solver\0";

  if (ksp_name == nullptr)
    ksp_name = NONAME_SOLVER;

  // Print message
  std::stringstream buff;
  buff << ksp_name << " iteration " << std::setw(4) << n << " - Residual " << std::scientific
       << std::setprecision(7) << rnorm / rhs_norm << std::endl;

  log.Log() << buff.str();

  return 0;
}

void
CopyVecToSTLvector(Vec x, std::vector<double>& data, size_t N, bool resize_STL)
{
  if (resize_STL)
  {
    data.clear();
    data.assign(N, 0.0);
  }
  else
    OpenSnLogicalErrorIf(data.size() < N,
                         "data.size() < N, " + std::to_string(data.size()) + " < " +
                           std::to_string(N));

  const double* x_ref = nullptr;
  VecGetArrayRead(x, &x_ref);

  std::copy(x_ref, x_ref + N, data.begin());

  VecRestoreArrayRead(x, &x_ref);
}

void
CopyVecToSTLvectorWithGhosts(Vec x, std::vector<double>& data, size_t N, bool resize_STL)
{
  if (resize_STL)
  {
    data.clear();
    data.assign(N, 0.0);
  }
  else
    OpenSnLogicalErrorIf(data.size() != N,
                         "data.size() != N, " + std::to_string(data.size()) + " < " +
                           std::to_string(N));

  auto info = GetGhostVectorLocalViewRead(x);
  const double* x_ref = info.x_localized_raw;

  std::copy(x_ref, x_ref + N, data.begin());

  RestoreGhostVectorLocalViewRead(x, info);
}

void
CopySTLvectorToVec(const std::vector<double>& data, Vec x, size_t N)
{
  double* x_ref = nullptr;
  VecGetArray(x, &x_ref);

  std::copy(data.begin(), data.end(), x_ref);

  VecRestoreArray(x, &x_ref);
}

void
CopyParallelVectorToVec(const ParallelVector& y, Vec x)
{
  const double* y_data = y.GetData();
  double* x_data = nullptr;
  VecGetArray(x, &x_data);
  std::copy(y_data, y_data + y.GetLocalSize(), x_data);
  VecRestoreArray(x, &x_data);
}

void
CopyGlobalVecToSTLvector(Vec x,
                         const std::vector<int64_t>& global_indices,
                         std::vector<double>& data)
{
  // Populating local indices
  auto N = static_cast<PetscInt>(global_indices.size());
  std::vector<int64_t> local_indices(N, 0);
  for (PetscInt counter = 0; counter < N; ++counter)
    local_indices[counter] = counter;

  // Creating PETSc vector
  Vec local_vec = nullptr;
  VecCreateSeq(PETSC_COMM_SELF, N + 1, &local_vec);
  VecSet(local_vec, 0.0);

  // Create and transfer index sets
  IS global_set = nullptr;
  IS local_set = nullptr;
  ISCreateGeneral(PETSC_COMM_SELF, N, global_indices.data(), PETSC_COPY_VALUES, &global_set);
  ISCreateGeneral(PETSC_COMM_SELF, N, local_indices.data(), PETSC_COPY_VALUES, &local_set);
  VecScatter scat = nullptr;
  VecScatterCreate(x, global_set, local_vec, local_set, &scat);
  VecScatterBegin(scat, x, local_vec, INSERT_VALUES, SCATTER_FORWARD);
  VecScatterEnd(scat, x, local_vec, INSERT_VALUES, SCATTER_FORWARD);

  // Copy to STL
  data.clear();
  data.resize(N, 0.0);
  const double* x_ref = nullptr;
  VecGetArrayRead(local_vec, &x_ref);

  std::copy(x_ref, x_ref + N, data.begin());

  VecRestoreArrayRead(x, &x_ref);

  // Cleanup
  ISDestroy(&global_set);
  ISDestroy(&local_set);

  VecDestroy(&local_vec);
}

void
CommunicateGhostEntries(Vec x)
{
  VecGhostUpdateBegin(x, INSERT_VALUES, SCATTER_FORWARD);
  VecGhostUpdateEnd(x, INSERT_VALUES, SCATTER_FORWARD);
}

GhostVecLocalRaw
GetGhostVectorLocalViewRead(Vec x)
{
  GhostVecLocalRaw local_data;

  VecGhostGetLocalForm(x, &local_data.x_localized);

  local_data.x_localized_raw = nullptr;
  VecGetArrayRead(local_data.x_localized, &local_data.x_localized_raw);

  return local_data;
}

void
RestoreGhostVectorLocalViewRead(Vec x, GhostVecLocalRaw& local_data)
{
  VecRestoreArrayRead(local_data.x_localized, &local_data.x_localized_raw);
  VecGhostRestoreLocalForm(x, &local_data.x_localized);
}

std::string
GetPETScConvergedReasonstring(KSPConvergedReason reason)
{
  std::stringstream ostr;
  switch (reason)
  {
    case KSP_CONVERGED_RTOL_NORMAL:
      ostr << "KSP_CONVERGED_RTOL_NORMAL";
      break;
    case KSP_CONVERGED_ATOL_NORMAL:
      ostr << "KSP_CONVERGED_ATOL_NORMAL";
      break;
    case KSP_CONVERGED_RTOL:
      ostr << "KSP_CONVERGED_RTOL";
      break;
    case KSP_CONVERGED_ATOL:
      ostr << "KSP_CONVERGED_ATOL";
      break;
    case KSP_CONVERGED_ITS:
      ostr << "KSP_CONVERGED_ITS";
      break;
#if PETSC_VERSION_LT(3, 19, 0)
    case KSP_CONVERGED_CG_NEG_CURVE:
      ostr << "KSP_CONVERGED_CG_NEG_CURVE";
      break;
#else
    case KSP_CONVERGED_NEG_CURVE:
      ostr << "KSP_CONVERGED_NEG_CURVE";
      break;
#endif
#if PETSC_VERSION_LT(3, 19, 0)
    case KSP_CONVERGED_CG_CONSTRAINED:
      ostr << "KSP_CONVERGED_CG_CONSTRAINED";
      break;
#endif
    case KSP_CONVERGED_STEP_LENGTH:
      ostr << "KSP_CONVERGED_STEP_LENGTH";
      break;
    case KSP_CONVERGED_HAPPY_BREAKDOWN:
      ostr << "KSP_CONVERGED_HAPPY_BREAKDOWN";
      break;
      /* diverged */
    case KSP_DIVERGED_NULL:
      ostr << "KSP_DIVERGED_NULL";
      break;
    case KSP_DIVERGED_ITS:
      ostr << "KSP_DIVERGED_ITS";
      break;
    case KSP_DIVERGED_DTOL:
      ostr << "KSP_DIVERGED_DTOL";
      break;
    case KSP_DIVERGED_BREAKDOWN:
      ostr << "KSP_DIVERGED_BREAKDOWN";
      break;
    case KSP_DIVERGED_BREAKDOWN_BICG:
      ostr << "KSP_DIVERGED_BREAKDOWN_BICG";
      break;
    case KSP_DIVERGED_NONSYMMETRIC:
      ostr << "KSP_DIVERGED_NONSYMMETRIC";
      break;
    case KSP_DIVERGED_INDEFINITE_PC:
      ostr << "KSP_DIVERGED_INDEFINITE_PC";
      break;
    case KSP_DIVERGED_NANORINF:
      ostr << "KSP_DIVERGED_NANORINF";
      break;
    case KSP_DIVERGED_INDEFINITE_MAT:
      ostr << "KSP_DIVERGED_INDEFINITE_MAT";
      break;

    default:
      ostr << "Unknown convergence reason.";
  }

  return ostr.str();
}

} // namespace opensn
