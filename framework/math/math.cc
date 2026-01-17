// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "framework/math/math.h"
#include "framework/runtime.h"
#include <cassert>
#include <petscmat.h>
#include <petscksp.h>
#include <iostream>
#include <iomanip>

namespace opensn
{

double
Factorial(const int x)
{
  double factorial_value = 1.0;
  for (int i = 2; i <= x; ++i)
    factorial_value *= i;

  return factorial_value;
}

std::pair<double, double>
OmegaToPhiThetaSafe(const Vector3& omega)
{
  // Notes: asin maps [-1,+1] to [-pi/2,+pi/2]
  //        acos maps [-1,+1] to [0,pi]
  // This mapping requires some logic for determining the azimuthal angle.
  //
  const auto omega_hat = omega.Normalized();

  double mu = omega_hat.z;
  mu = std::min(mu, 1.0);
  mu = std::max(mu, -1.0);

  double theta = acos(mu);

  // Handling omega aligned to k_hat
  if (std::fabs(omega_hat.z) < 1.0e-16)
    return {0.0, theta};

  double cos_phi = omega_hat.x / sin(theta);
  cos_phi = std::min(cos_phi, 1.0);
  cos_phi = std::max(cos_phi, -1.0);

  // Computing varphi for NE and NW quadrant
  if (omega_hat.y >= 0.0)
    return {acos(cos_phi), theta};
  // Computing varphi for SE and SW quadrant
  else
    return {2.0 * M_PI - acos(cos_phi), theta};
}

std::vector<std::vector<double>>
Transpose(const std::vector<std::vector<double>>& matrix)
{
  // Check if matrix is empty
  if (matrix.empty())
  {
    throw std::runtime_error("Cannot transpose empty matrix");
  }

  const size_t m = matrix.size();    // number of rows
  const size_t n = matrix[0].size(); // number of columns

  // Verify consistent row dimensions
  for (size_t i = 1; i < m; ++i)
  {
    if (matrix[i].size() != n)
    {
      throw std::runtime_error("Matrix must have consistent row dimensions");
    }
  }

  // Create transposed matrix with dimensions n x m
  std::vector<std::vector<double>> transposed(n, std::vector<double>(m));

  // Perform transpose: transposed[j][i] = matrix[i][j]
  for (size_t i = 0; i < m; ++i)
  {
    for (size_t j = 0; j < n; ++j)
    {
      transposed[j][i] = matrix[i][j];
    }
  }

  return transposed;
}

std::vector<std::vector<double>>
InvertMatrix(const std::vector<std::vector<double>>& matrix)
{
  const PetscInt n = static_cast<PetscInt>(matrix.size());

  // Check if matrix is square
  if (n == 0)
  {
    throw std::runtime_error("Matrix must not be empty for inversion");
  }

  for (PetscInt i = 0; i < n; ++i)
  {
    if (static_cast<PetscInt>(matrix[i].size()) != n)
    {
      throw std::runtime_error("Matrix must be square for inversion");
    }
  }

  Mat A, B, X;
  MatFactorInfo info;

  // Create dense matrix A directly
  MatCreateSeqDense(PETSC_COMM_SELF, n, n, NULL, &A);

  // Fill matrix A using array access (column-major storage)
  PetscScalar* aArray;
  MatDenseGetArrayWrite(A, &aArray);

  for (PetscInt j = 0; j < n; ++j) // columns
  {
    for (PetscInt i = 0; i < n; ++i) // rows
    {
      aArray[j * n + i] = matrix[i][j]; // Column-major: col j, row i
    }
  }

  MatDenseRestoreArrayWrite(A, &aArray);
  MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY);

  // Create identity matrix B (right-hand side)
  MatCreateSeqDense(PETSC_COMM_SELF, n, n, NULL, &B);

  // Fill B with identity matrix using array access
  PetscScalar* bArray;
  MatDenseGetArrayWrite(B, &bArray);

  // Initialize to zero first
  for (PetscInt idx = 0; idx < n * n; ++idx)
  {
    bArray[idx] = 0.0;
  }
  // Set diagonal to 1
  for (PetscInt i = 0; i < n; ++i)
  {
    bArray[i * n + i] = 1.0; // Column-major: diagonal elements
  }

  MatDenseRestoreArrayWrite(B, &bArray);
  MatAssemblyBegin(B, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(B, MAT_FINAL_ASSEMBLY);

  // Create solution matrix X
  MatDuplicate(B, MAT_DO_NOT_COPY_VALUES, &X);

  // Perform LU factorization directly
  IS rowPerm, colPerm;
  Mat F; // Factored matrix

  MatGetOrdering(A, MATORDERINGND, &rowPerm, &colPerm);
  MatFactorInfoInitialize(&info);
  MatGetFactor(A, MATSOLVERPETSC, MAT_FACTOR_LU, &F);
  MatLUFactorSymbolic(F, A, rowPerm, colPerm, &info);

  // This is the critical operation that can fail for singular matrices
  PetscErrorCode ierr = MatLUFactorNumeric(F, A, &info);
  if (ierr != 0)
  {
    // Clean up before throwing
    MatDestroy(&A);
    MatDestroy(&B);
    MatDestroy(&X);
    MatDestroy(&F);
    ISDestroy(&rowPerm);
    ISDestroy(&colPerm);
    throw std::runtime_error("LU factorization failed - matrix may be singular");
  }

  // Check if factorization succeeded (check for zero pivot)
  PetscBool zeropivot;
  MatFactorGetError(F, (MatFactorError*)&zeropivot);

  if (zeropivot)
  {
    // Clean up before throwing
    MatDestroy(&A);
    MatDestroy(&B);
    MatDestroy(&X);
    MatDestroy(&F);
    ISDestroy(&rowPerm);
    ISDestroy(&colPerm);
    throw std::runtime_error(
      "Matrix is singular or nearly singular (zero pivot in LU factorization)");
  }

  // Solve AX = B (where B is identity, so X will be A^-1)
  MatMatSolve(F, B, X);

  // Extract solution back to std::vector format
  std::vector<std::vector<double>> inverse(n, std::vector<double>(n));

  // Get array from PETSc dense matrix (read-only)
  const PetscScalar* xArray;
  MatDenseGetArrayRead(X, &xArray);

  // Copy data (PETSc uses column-major storage for dense matrices)
  for (PetscInt i = 0; i < n; ++i)
  {
    for (PetscInt j = 0; j < n; ++j)
    {
      inverse[i][j] = PetscRealPart(xArray[j * n + i]); // Transpose due to column-major
    }
  }

  MatDenseRestoreArrayRead(X, &xArray);

  // Clean up PETSc objects
  MatDestroy(&A);
  MatDestroy(&B);
  MatDestroy(&X);
  MatDestroy(&F);
  ISDestroy(&rowPerm);
  ISDestroy(&colPerm);

  return inverse;
}

std::vector<std::vector<double>>
OrthogonalizeHouseholder(const std::vector<std::vector<double>>& matrix,
                         const std::vector<double>& weights)
{
  // Check an empty matrix
  if (matrix.empty())
  {
    throw std::runtime_error("Cannot orthogonalize empty matrix");
  }

  const size_t m = matrix.size();
  const size_t n = matrix[0].size();

  // Verify rectangular matrix
  for (size_t i = 1; i < m; ++i)
  {
    if (matrix[i].size() != n)
    {
      throw std::runtime_error("Matrix must have consistent column dimensions");
    }
  }

  // Orthogonalize columns using Modified Gram-Schmidt with weighted inner product
  // This keeps each column as close to its original form as possible while ensuring orthogonality

  // Create working copy of the matrix - we will orthogonalize its columns in place
  std::vector<std::vector<double>> Q = matrix;

  // Orthogonalize columns using Modified Gram-Schmidt
  for (size_t j = 0; j < n; ++j)
  {
    // Orthogonalize column j against all previous columns
    for (size_t i = 0; i < j; ++i)
    {
      // Compute weighted inner product <Q[:,i], Q[:,j]>
      double dot_product = 0.0;
      for (size_t row = 0; row < m; ++row)
      {
        dot_product += weights[row] * Q[row][i] * Q[row][j];
      }

      // Subtract projection: Q[:,j] -= dot_product * Q[:,i]
      for (size_t row = 0; row < m; ++row)
      {
        Q[row][j] -= dot_product * Q[row][i];
      }
    }

    // Normalize column j under weighted inner product
    double norm_squared = 0.0;
    for (size_t row = 0; row < m; ++row)
    {
      norm_squared += weights[row] * Q[row][j] * Q[row][j];
    }
    double norm = std::sqrt(norm_squared);

    if (norm > 1e-14)
    {
      for (size_t row = 0; row < m; ++row)
      {
        Q[row][j] /= norm;
      }
    }
  }
  return Q;
}

void
PrintVector(const std::vector<double>& x)
{
  for (const auto& xi : x)
    std::cout << xi << ' ';
  std::cout << std::endl;
}

void
Scale(std::vector<double>& x, const double& val)
{
  for (double& xi : x)
    xi *= val;
}

void
Set(std::vector<double>& x, const double& val)
{
  for (double& xi : x)
    xi = val;
}

double
Dot(const std::vector<double>& x, const std::vector<double>& y)
{
  // Error Checks
  assert(!x.empty());
  assert(!y.empty());
  assert(x.size() == y.size());
  // Local Variables
  size_t n = x.size();
  double val = 0.0;

  for (size_t i = 0; i != n; ++i)
    val += x[i] * y[i];

  return val;
}

std::vector<double>
Mult(const std::vector<double>& x, const double& val)
{
  size_t n = x.size();
  std::vector<double> y(n);

  for (size_t i = 0; i != n; ++i)
    y[i] = val * x[i];

  return y;
}

double
L1Norm(const std::vector<double>& x)
{
  // Local Variables
  size_t n = x.size();
  double val = 0.0;

  for (size_t i = 0; i != n; ++i)
    val += std::fabs(x[i]);

  return val;
}

double
L2Norm(const std::vector<double>& x)
{
  // Local Variables
  size_t n = x.size();
  double val = 0.0;

  for (size_t i = 0; i != n; ++i)
    val += x[i] * x[i];

  return std::sqrt(val);
}

double
LInfNorm(const std::vector<double>& x)
{
  // Local Variables
  size_t n = x.size();
  double val = 0.0;

  for (size_t i = 0; i != n; ++i)
    val += std::max(std::fabs(x[i]), val);

  return val;
}

double
LpNorm(const std::vector<double>& x, const double& p)
{
  // Local Variables
  size_t n = x.size();
  double val = 0.0;

  for (size_t i = 0; i != n; ++i)
    val += std::pow(std::fabs(x[i]), p);

  return std::pow(val, 1.0 / p);
}

std::vector<double>
operator+(const std::vector<double>& a, const std::vector<double>& b)
{
  assert(a.size() == b.size());
  std::vector<double> result(a.size(), 0.0);

  for (size_t i = 0; i < a.size(); ++i)
    result[i] = a[i] + b[i];

  return result;
}

std::vector<double>
operator-(const std::vector<double>& a, const std::vector<double>& b)
{
  assert(a.size() == b.size());
  std::vector<double> result(a.size(), 0.0);

  for (size_t i = 0; i < a.size(); ++i)
    result[i] = a[i] - b[i];

  return result;
}

double
ComputeL2Change(std::vector<double>& x, std::vector<double>& y)
{
  assert(x.size() == y.size());

  double norm = 0.0;
  for (auto i = 0; i < x.size(); ++i)
  {
    double val = x[i] - y[i];
    norm += val * val;
  }

  double global_norm = 0.0;
  mpi_comm.all_reduce<double>(norm, global_norm, mpi::op::sum<double>());

  return std::sqrt(global_norm);
}

double
ComputePointwiseChange(std::vector<double>& x, std::vector<double>& y)
{
  assert(x.size() == y.size());

  double pw_change = 0.0;
  for (auto i = 0; i < x.size(); ++i)
  {
    double max = std::max(x[i], y[i]);
    double delta = std::fabs(x[i] - y[i]);
    if (max >= std::numeric_limits<double>::min())
      pw_change = std::max(delta / max, pw_change);
    else
      pw_change = std::max(delta, pw_change);
  }

  double global_pw_change = 0.0;
  mpi_comm.all_reduce<double>(pw_change, global_pw_change, mpi::op::max<double>());

  return global_pw_change;
}

} // namespace opensn
