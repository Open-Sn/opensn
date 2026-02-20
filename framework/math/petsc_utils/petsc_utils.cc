// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "framework/math/petsc_utils/petsc_utils.h"
#include "framework/data_types/parallel_vector/parallel_vector.h"
#include "framework/runtime.h"
#include "framework/logging/log.h"
#include <iomanip>
#include <sstream>
#include <limits>

namespace opensn
{
namespace
{
PetscInt
ToPetscInt(const int64_t value, const char* name)
{
  OpenSnLogicalErrorIf(value < static_cast<int64_t>(std::numeric_limits<PetscInt>::min()) or
                         value > static_cast<int64_t>(std::numeric_limits<PetscInt>::max()),
                       std::string(name) + " is out of PetscInt range");
  return static_cast<PetscInt>(value);
}
} // namespace

[[noreturn]] void
ThrowPETScError(int ierr, const char* expr, const char* file, int line)
{
  const char* ierr_desc = nullptr;
  char* ierr_desc_specific = nullptr;
  PetscErrorMessage(static_cast<PetscErrorCode>(ierr), &ierr_desc, &ierr_desc_specific);

  std::stringstream ss;
  ss << "PETSc call failed at " << file << ":" << line << " in expression \"" << expr
     << "\" with error code " << ierr;
  if (ierr_desc != nullptr)
    ss << ". " << ierr_desc;
  if (ierr_desc_specific != nullptr)
    ss << " (" << ierr_desc_specific << ")";

  throw std::runtime_error(ss.str());
}

Vec
CreateVector(int64_t local_size, int64_t global_size)
{
  Vec x = nullptr;
  OpenSnPETScCall(VecCreate(opensn::mpi_comm, &x));
  OpenSnPETScCall(VecSetType(x, VECMPI));
  OpenSnPETScCall(
    VecSetSizes(x, ToPetscInt(local_size, "local_size"), ToPetscInt(global_size, "global_size")));
  OpenSnPETScCall(VecSetOption(x, VEC_IGNORE_NEGATIVE_INDICES, PETSC_TRUE));

  return x;
}

void
CreateVector(Vec& x, int64_t local_size, int64_t global_size)
{
  OpenSnPETScCall(VecCreate(opensn::mpi_comm, &x));
  OpenSnPETScCall(VecSetType(x, VECMPI));
  OpenSnPETScCall(
    VecSetSizes(x, ToPetscInt(local_size, "local_size"), ToPetscInt(global_size, "global_size")));
  OpenSnPETScCall(VecSetOption(x, VEC_IGNORE_NEGATIVE_INDICES, PETSC_TRUE));
}

Vec
CreateVectorWithGhosts(int64_t local_size,
                       int64_t global_size,
                       int64_t nghosts,
                       const std::vector<int64_t>& ghost_indices)
{
  std::vector<PetscInt> ghost_indices_petsc;
  ghost_indices_petsc.reserve(ghost_indices.size());
  for (const auto gid : ghost_indices)
    ghost_indices_petsc.push_back(ToPetscInt(gid, "ghost_index"));

  Vec x = nullptr;
  OpenSnPETScCall(
    VecCreateGhost(opensn::mpi_comm,
                   ToPetscInt(local_size, "local_size"),
                   ToPetscInt(global_size, "global_size"),
                   ToPetscInt(nghosts, "nghosts"),
                   (ghost_indices_petsc.empty()) ? nullptr : ghost_indices_petsc.data(),
                   &x));

  OpenSnPETScCall(VecSetOption(x, VEC_IGNORE_NEGATIVE_INDICES, PETSC_TRUE));

  return x;
}

Mat
CreateSquareMatrix(int64_t local_size, int64_t global_size)
{
  Mat A = nullptr;
  OpenSnPETScCall(MatCreate(opensn::mpi_comm, &A));
  OpenSnPETScCall(MatSetType(A, MATMPIAIJ));
  const auto local = ToPetscInt(local_size, "local_size");
  const auto global = ToPetscInt(global_size, "global_size");
  OpenSnPETScCall(MatSetSizes(A, local, local, global, global));

  OpenSnPETScCall(MatMPIAIJSetPreallocation(A, 1, nullptr, 0, nullptr));
  OpenSnPETScCall(MatSetOption(A, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE));
  OpenSnPETScCall(MatSetOption(A, MAT_IGNORE_ZERO_ENTRIES, PETSC_TRUE));

  return A;
}

void
CreateSquareMatrix(Mat& A, int64_t local_size, int64_t global_size)
{
  OpenSnPETScCall(MatCreate(opensn::mpi_comm, &A));
  OpenSnPETScCall(MatSetType(A, MATMPIAIJ));
  const auto local = ToPetscInt(local_size, "local_size");
  const auto global = ToPetscInt(global_size, "global_size");
  OpenSnPETScCall(MatSetSizes(A, local, local, global, global));

  OpenSnPETScCall(MatMPIAIJSetPreallocation(A, 1, nullptr, 0, nullptr));
  OpenSnPETScCall(MatSetOption(A, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE));
  OpenSnPETScCall(MatSetOption(A, MAT_IGNORE_ZERO_ENTRIES, PETSC_TRUE));
}

void
InitMatrixSparsity(Mat& A,
                   const std::vector<int64_t>& nodal_nnz_in_diag,
                   const std::vector<int64_t>& nodal_nnz_off_diag)
{
  std::vector<PetscInt> nodal_nnz_in_diag_petsc;
  std::vector<PetscInt> nodal_nnz_off_diag_petsc;
  nodal_nnz_in_diag_petsc.reserve(nodal_nnz_in_diag.size());
  nodal_nnz_off_diag_petsc.reserve(nodal_nnz_off_diag.size());
  for (const auto nnz : nodal_nnz_in_diag)
    nodal_nnz_in_diag_petsc.push_back(ToPetscInt(nnz, "nodal_nnz_in_diag"));
  for (const auto nnz : nodal_nnz_off_diag)
    nodal_nnz_off_diag_petsc.push_back(ToPetscInt(nnz, "nodal_nnz_off_diag"));

  OpenSnPETScCall(MatMPIAIJSetPreallocation(
    A, 0, nodal_nnz_in_diag_petsc.data(), 0, nodal_nnz_off_diag_petsc.data()));
  OpenSnPETScCall(MatSetOption(A, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE));
  OpenSnPETScCall(MatSetOption(A, MAT_IGNORE_ZERO_ENTRIES, PETSC_TRUE));
  OpenSnPETScCall(MatSetUp(A));
}

void
InitMatrixSparsity(Mat& A, int64_t nodal_nnz_in_diag, int64_t nodal_nnz_off_diag)
{
  OpenSnPETScCall(MatMPIAIJSetPreallocation(A,
                                            ToPetscInt(nodal_nnz_in_diag, "nodal_nnz_in_diag"),
                                            nullptr,
                                            ToPetscInt(nodal_nnz_off_diag, "nodal_nnz_off_diag"),
                                            nullptr));
  OpenSnPETScCall(MatSetOption(A, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE));
  OpenSnPETScCall(MatSetOption(A, MAT_IGNORE_ZERO_ENTRIES, PETSC_TRUE));
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

  OpenSnPETScCall(KSPCreate(opensn::mpi_comm, &setup.ksp));
  OpenSnPETScCall(KSPSetOperators(setup.ksp, matrix, matrix));
  OpenSnPETScCall(KSPSetType(setup.ksp, solver_type.c_str()));

  OpenSnPETScCall(KSPSetOptionsPrefix(setup.ksp, solver_name.c_str()));

  OpenSnPETScCall(KSPGetPC(setup.ksp, &setup.pc));
  OpenSnPETScCall(PCSetType(setup.pc, preconditioner_type.c_str()));

  OpenSnPETScCall(KSPSetTolerances(
    setup.ksp, rel_tol, abs_tol, 1.0e50, ToPetscInt(maximum_iterations, "maximum_iterations")));
  OpenSnPETScCall(KSPSetInitialGuessNonzero(setup.ksp, PETSC_TRUE));

  OpenSnPETScCall(KSPSetFromOptions(setup.ksp));

  OpenSnPETScCall(KSPMonitorSet(setup.ksp, &KSPMonitorRelativeToRHS, nullptr, nullptr));

  return setup;
}

PetscErrorCode
KSPMonitorRelativeToRHS(KSP ksp, PetscInt n, PetscReal rnorm, void* /* context */)
{
  Vec Rhs = nullptr;
  PetscErrorCode ierr = KSPGetRhs(ksp, &Rhs);
  if (ierr != PETSC_SUCCESS)
    return ierr;
  double rhs_norm = 0.0;
  ierr = VecNorm(Rhs, NORM_2, &rhs_norm);
  if (ierr != PETSC_SUCCESS)
    return ierr;
  if (rhs_norm < 1.0e-12)
    rhs_norm = 1.0;

  // Get solver name
  const char* ksp_name = nullptr;
  ierr = KSPGetOptionsPrefix(ksp, &ksp_name);
  if (ierr != PETSC_SUCCESS)
    return ierr;

  // Default to this if ksp_name is NULL
  const char NONAME_SOLVER[] = "NoName-Solver\0";

  if (ksp_name == nullptr)
    ksp_name = NONAME_SOLVER;

  // Print message
  std::stringstream buff;
  buff << ksp_name << " iteration " << std::setw(4) << n << " - Residual " << std::scientific
       << std::setprecision(7) << rnorm / rhs_norm << std::endl;

  log.Log() << buff.str();

  return PETSC_SUCCESS;
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
  OpenSnPETScCall(VecGetArrayRead(x, &x_ref));

  std::copy(x_ref, x_ref + N, data.begin());

  OpenSnPETScCall(VecRestoreArrayRead(x, &x_ref));
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
  OpenSnPETScCall(VecGetArray(x, &x_ref));

  std::copy(data.begin(), data.end(), x_ref);

  OpenSnPETScCall(VecRestoreArray(x, &x_ref));
}

void
CopyParallelVectorToVec(const ParallelVector& y, Vec x)
{
  const double* y_data = y.GetData();
  double* x_data = nullptr;
  OpenSnPETScCall(VecGetArray(x, &x_data));
  std::copy(y_data, y_data + y.GetLocalSize(), x_data);
  OpenSnPETScCall(VecRestoreArray(x, &x_data));
}

void
CopyGlobalVecToSTLvector(Vec x,
                         const std::vector<int64_t>& global_indices,
                         std::vector<double>& data)
{
  std::vector<PetscInt> global_indices_petsc;
  global_indices_petsc.reserve(global_indices.size());
  for (const auto idx : global_indices)
    global_indices_petsc.push_back(ToPetscInt(idx, "global_index"));

  // Populating local indices
  const auto N = static_cast<PetscInt>(global_indices_petsc.size());
  std::vector<PetscInt> local_indices(N, 0);
  for (PetscInt counter = 0; counter < N; ++counter)
    local_indices[counter] = counter;

  // Creating PETSc vector
  Vec local_vec = nullptr;
  OpenSnPETScCall(VecCreateSeq(PETSC_COMM_SELF, N + 1, &local_vec));
  OpenSnPETScCall(VecSet(local_vec, 0.0));

  // Create and transfer index sets
  IS global_set = nullptr;
  IS local_set = nullptr;
  OpenSnPETScCall(ISCreateGeneral(
    PETSC_COMM_SELF, N, global_indices_petsc.data(), PETSC_COPY_VALUES, &global_set));
  OpenSnPETScCall(
    ISCreateGeneral(PETSC_COMM_SELF, N, local_indices.data(), PETSC_COPY_VALUES, &local_set));
  VecScatter scat = nullptr;
  OpenSnPETScCall(VecScatterCreate(x, global_set, local_vec, local_set, &scat));
  OpenSnPETScCall(VecScatterBegin(scat, x, local_vec, INSERT_VALUES, SCATTER_FORWARD));
  OpenSnPETScCall(VecScatterEnd(scat, x, local_vec, INSERT_VALUES, SCATTER_FORWARD));

  // Copy to STL
  data.clear();
  data.resize(N, 0.0);
  const double* x_ref = nullptr;
  OpenSnPETScCall(VecGetArrayRead(local_vec, &x_ref));

  std::copy(x_ref, x_ref + N, data.begin());

  OpenSnPETScCall(VecRestoreArrayRead(x, &x_ref));

  // Cleanup
  OpenSnPETScCall(ISDestroy(&global_set));
  OpenSnPETScCall(ISDestroy(&local_set));

  OpenSnPETScCall(VecDestroy(&local_vec));
}

void
CommunicateGhostEntries(Vec x)
{
  OpenSnPETScCall(VecGhostUpdateBegin(x, INSERT_VALUES, SCATTER_FORWARD));
  OpenSnPETScCall(VecGhostUpdateEnd(x, INSERT_VALUES, SCATTER_FORWARD));
}

GhostVecLocalRaw
GetGhostVectorLocalViewRead(Vec x)
{
  GhostVecLocalRaw local_data;

  OpenSnPETScCall(VecGhostGetLocalForm(x, &local_data.x_localized));

  local_data.x_localized_raw = nullptr;
  OpenSnPETScCall(VecGetArrayRead(local_data.x_localized, &local_data.x_localized_raw));

  return local_data;
}

void
RestoreGhostVectorLocalViewRead(Vec x, GhostVecLocalRaw& local_data)
{
  OpenSnPETScCall(VecRestoreArrayRead(local_data.x_localized, &local_data.x_localized_raw));
  OpenSnPETScCall(VecGhostRestoreLocalForm(x, &local_data.x_localized));
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
