#include "framework/math/petsc_utils/petsc_utils.h"
#include "framework/math/parallel_vector/parallel_vector.h"
#include "framework/runtime.h"
#include "framework/logging/log.h"
#include <iomanip>

namespace opensn
{

Vec
CreateVector(int64_t local_size, int64_t global_size)
{
  Vec x;
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
                       const std::vector<int64_t>& ghost_indices)
{
  Vec x;
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
  Mat A;
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
CreateCommonKrylovSolverSetup(Mat ref_matrix,
                              const std::string& in_solver_name,
                              const std::string& in_solver_type,
                              const std::string& in_preconditioner_type,
                              double in_relative_residual_tolerance,
                              int64_t in_maximum_iterations)
{
  PETScSolverSetup setup;

  KSPCreate(opensn::mpi_comm, &setup.ksp);
  KSPSetOperators(setup.ksp, ref_matrix, ref_matrix);
  KSPSetType(setup.ksp, in_solver_type.c_str());

  KSPSetOptionsPrefix(setup.ksp, in_solver_name.c_str());

  KSPGetPC(setup.ksp, &setup.pc);
  PCSetType(setup.pc, in_preconditioner_type.c_str());

  KSPSetTolerances(
    setup.ksp, 1.e-50, in_relative_residual_tolerance, 1.0e50, in_maximum_iterations);
  KSPSetInitialGuessNonzero(setup.ksp, PETSC_TRUE);

  KSPSetFromOptions(setup.ksp);

  KSPMonitorSet(setup.ksp, &KSPMonitorRelativeToRHS, nullptr, nullptr);

  return setup;
}

PetscErrorCode
KSPMonitorRelativeToRHS(KSP ksp, PetscInt n, PetscReal rnorm, void*)
{
  Vec Rhs;
  KSPGetRhs(ksp, &Rhs);
  double rhs_norm;
  VecNorm(Rhs, NORM_2, &rhs_norm);
  if (rhs_norm < 1.0e-12)
    rhs_norm = 1.0;

  // Get solver name
  const char* ksp_name;
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

PetscErrorCode
KSPMonitorStraight(KSP ksp, PetscInt n, PetscReal rnorm, void*)
{
  // Get solver name
  const char* ksp_name;
  KSPGetOptionsPrefix(ksp, &ksp_name);

  // Default to this if ksp_name is NULL
  const char NONAME_SOLVER[] = "NoName-Solver\0";

  if (ksp_name == nullptr)
    ksp_name = NONAME_SOLVER;

  // Print message
  std::stringstream buff;
  buff << ksp_name << " iteration " << std::setw(4) << n << " - Residual " << std::scientific
       << std::setprecision(7) << rnorm << std::endl;

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
    ChiLogicalErrorIf(data.size() < N,
                      "data.size() < N, " + std::to_string(data.size()) + " < " +
                        std::to_string(N));

  const double* x_ref;
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
    ChiLogicalErrorIf(data.size() != N,
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
  double* x_ref;
  VecGetArray(x, &x_ref);

  std::copy(data.begin(), data.end(), x_ref);

  VecRestoreArray(x, &x_ref);
}

void
CopyParallelVectorToVec(const ParallelVector& y, Vec x)
{
  const double* y_data = y.Data();
  double* x_data;
  VecGetArray(x, &x_data);
  std::copy(y_data, y_data + y.LocalSize(), x_data);
  VecRestoreArray(x, &x_data);
}

void
CopyGlobalVecToSTLvector(Vec x,
                         const std::vector<int64_t>& global_indices,
                         std::vector<double>& data)
{
  // Populating local indices
  size_t N = global_indices.size();
  std::vector<int64_t> local_indices(N, 0);
  size_t counter = 0;
  for (auto val : global_indices)
  {
    local_indices[counter] = counter;
    ++counter;
  }

  // Creating PETSc vector
  Vec local_vec;
  VecCreateSeq(PETSC_COMM_SELF, global_indices.size() + 1, &local_vec);
  VecSet(local_vec, 0.0);

  // Create and transfer index sets
  IS global_set;
  IS local_set;
  ISCreateGeneral(PETSC_COMM_SELF, N, global_indices.data(), PETSC_COPY_VALUES, &global_set);
  ISCreateGeneral(PETSC_COMM_SELF, N, local_indices.data(), PETSC_COPY_VALUES, &local_set);
  VecScatter scat;
  VecScatterCreate(x, global_set, local_vec, local_set, &scat);
  VecScatterBegin(scat, x, local_vec, INSERT_VALUES, SCATTER_FORWARD);
  VecScatterEnd(scat, x, local_vec, INSERT_VALUES, SCATTER_FORWARD);

  // Copy to STL
  data.clear();
  data.resize(N, 0.0);
  const double* x_ref;
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
  Vec x_localized;
  VecGhostGetLocalForm(x, &x_localized);
  const double* x_localized_raw;

  VecGetArrayRead(x_localized, &x_localized_raw);

  GhostVecLocalRaw local_data;
  local_data.x_localized = x_localized;
  local_data.x_localized_raw = (double*)x_localized_raw;

  return local_data;
}

void
RestoreGhostVectorLocalViewRead(Vec x, GhostVecLocalRaw& local_data)
{
  VecRestoreArrayRead(local_data.x_localized, (const double**)&local_data.x_localized_raw);
  VecGhostRestoreLocalForm(x, &local_data.x_localized);
}

} // namespace opensn
