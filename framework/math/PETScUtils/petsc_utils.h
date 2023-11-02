#ifndef CHI_MATH_PETSC_UTILS_H
#define CHI_MATH_PETSC_UTILS_H

#include <petscksp.h>
#include <vector>

namespace chi_math
{
class ParallelVector;
}

namespace chi_math::PETScUtils
{
/**Generalized solver structure.*/
struct PETScSolverSetup
{
  KSP ksp;
  PC pc;

  std::string in_solver_name = "KSPSolver";

  std::string solver_type = KSPGMRES;
  std::string preconditioner_type = PCNONE;

  double relative_residual_tol = 1.0e-6;
  int maximum_iterations = 100;
};

// 01
Vec CreateVector(int64_t local_size, int64_t global_size);
void CreateVector(Vec& x, int64_t local_size, int64_t global_size);

Vec CreateVectorWithGhosts(int64_t local_size,
                           int64_t global_size,
                           int64_t nghosts,
                           const std::vector<int64_t>& ghost_indices);

// 02
Mat CreateSquareMatrix(int64_t local_size, int64_t global_size);
void CreateSquareMatrix(Mat& A, int64_t local_size, int64_t global_size);
void InitMatrixSparsity(Mat& A,
                        const std::vector<int64_t>& nodal_nnz_in_diag,
                        const std::vector<int64_t>& nodal_nnz_off_diag);
void InitMatrixSparsity(Mat& A, int64_t nodal_nnz_in_diag, int64_t nodal_nnz_off_diag);

// 03
PETScSolverSetup CreateCommonKrylovSolverSetup(Mat ref_matrix,
                                               const std::string& in_solver_name = "KSPSolver",
                                               const std::string& in_solver_type = KSPGMRES,
                                               const std::string& in_preconditioner_type = PCNONE,
                                               double in_relative_residual_tolerance = 1.0e-6,
                                               int64_t in_maximum_iterations = 100);

PetscErrorCode RelativeResidualConvergenceTest(
  KSP ksp, PetscInt n, PetscReal rnorm, KSPConvergedReason* convergedReason, void* monitordestroy);

PetscErrorCode KSPMonitorRelativeToRHS(KSP ksp, PetscInt n, PetscReal rnorm, void*);
PetscErrorCode KSPMonitorStraight(KSP ksp, PetscInt n, PetscReal rnorm, void*);

// 04
void CopyVecToSTLvector(Vec x, std::vector<double>& data, size_t N, bool resize_STL = true);
void
CopyVecToSTLvectorWithGhosts(Vec x, std::vector<double>& data, size_t N, bool resize_STL = true);

void CopySTLvectorToVec(const std::vector<double>& data, Vec x, size_t N);
void CopyParallelVectorToVec(const ParallelVector& y, Vec x);

void CopyGlobalVecToSTLvector(Vec x,
                              const std::vector<int64_t>& global_indices,
                              std::vector<double>& data);

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

GhostVecLocalRaw GetGhostVectorLocalViewRead(Vec x);
void RestoreGhostVectorLocalViewRead(Vec x, GhostVecLocalRaw& local_data);

} // namespace chi_math::PETScUtils

#endif