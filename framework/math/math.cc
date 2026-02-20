// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "framework/math/math.h"
#include "framework/runtime.h"
#include "framework/logging/log.h"
#include <cassert>
#include <petscmat.h>
#include <petscksp.h>
#include <petscblaslapack.h>
#include <iostream>
#include <iomanip>
#include <limits>

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
  const auto n = static_cast<PetscInt>(matrix.size());

  // Check if matrix is square
  if (n == 0)
  {
    throw std::runtime_error("Matrix must not be empty for inversion");
  }

  for (PetscInt i = 0; i < n; ++i)
  {
    if (std::cmp_not_equal(matrix[i].size(), n))
    {
      throw std::runtime_error("Matrix must be square for inversion");
    }
  }

  // Use SVD-based pseudo-inverse for numerical stability
  // This is more robust than LU factorization for ill-conditioned matrices

  // Create working arrays for LAPACK SVD (dgesvd)
  // A = U * S * V^T, so A^(-1) = V * S^(-1) * U^T
  std::vector<double> A_col_major(n * n);
  std::vector<double> U(n * n);
  std::vector<double> VT(n * n);
  std::vector<double> S(n);

  // Copy matrix to column-major format for LAPACK
  for (PetscInt j = 0; j < n; ++j)
  {
    for (PetscInt i = 0; i < n; ++i)
    {
      A_col_major[j * n + i] = matrix[i][j];
    }
  }

  // Call LAPACK dgesvd for SVD decomposition
  // jobu = 'A' means compute all columns of U
  // jobvt = 'A' means compute all rows of V^T
  char jobu = 'A';
  char jobvt = 'A';
  auto n_blas = static_cast<PetscBLASInt>(n);
  PetscBLASInt lwork = std::max(PetscBLASInt(1), 5 * n_blas);
  std::vector<double> work(lwork);
  PetscBLASInt info_svd = 0;

  // LAPACK SVD: dgesvd
  LAPACKgesvd_(&jobu,
               &jobvt,
               &n_blas,
               &n_blas,
               A_col_major.data(),
               &n_blas,
               S.data(),
               U.data(),
               &n_blas,
               VT.data(),
               &n_blas,
               work.data(),
               &lwork,
               &info_svd);

  if (info_svd != 0)
  {
    throw std::runtime_error("SVD decomposition failed with LAPACK error code: " +
                             std::to_string(info_svd));
  }

  // Compute condition number and determine tolerance for singular values
  double sigma_max = S[0];
  double sigma_min = S[n - 1];

  // Use machine epsilon scaled by matrix size and largest singular value as tolerance
  // This is a standard approach for pseudo-inverse computation
  const double machine_eps = std::numeric_limits<double>::epsilon();
  const double tol = static_cast<double>(n) * machine_eps * sigma_max;

  // Log condition number for debugging ill-conditioned matrices
  double condition_number = (sigma_min > tol) ? (sigma_max / sigma_min) : INFINITY;
  if (condition_number > 1.0e10)
  {
    log.Log0Warning() << "InvertMatrix: Matrix is ill-conditioned (condition number = "
                      << std::scientific << std::setprecision(3) << condition_number << ")";
  }

  // Compute pseudo-inverse: A^(-1) = V * S^(-1) * U^T
  // Since we have V^T from SVD, we need V = (V^T)^T
  // Result[i][j] = sum_k V[i][k] * (1/S[k]) * U[j][k]
  //              = sum_k VT[k][i] * (1/S[k]) * U[j][k]
  // In column-major: VT[k*n + i], U[k*n + j]

  std::vector<std::vector<double>> inverse(n, std::vector<double>(n, 0.0));

  for (PetscInt i = 0; i < n; ++i)
  {
    for (PetscInt j = 0; j < n; ++j)
    {
      double sum = 0.0;
      for (PetscInt k = 0; k < n; ++k)
      {
        // Only include contributions from singular values above tolerance
        if (S[k] > tol)
        {
          // V[i][k] = VT[k][i] (in row-major) = VT[i + k*n] (in col-major storage of VT)
          // U[j][k] = U[k][j]^T, but U is stored col-major, so U[j][k] = U[j + k*n]
          // Wait, let me reconsider the indexing carefully:
          // VT is stored column-major, VT(i,j) = VT[i + j*n]
          // U is stored column-major, U(i,j) = U[i + j*n]
          // We want: V(i,k) * U^T(k,j) = V(i,k) * U(j,k)
          // V = VT^T, so V(i,k) = VT(k,i) = VT[k + i*n]
          // U(j,k) = U[j + k*n]
          double v_ik = VT[k + i * n];
          double u_jk = U[j + k * n];
          sum += v_ik * (1.0 / S[k]) * u_jk;
        }
      }
      inverse[i][j] = sum;
    }
  }

  return inverse;
}

std::vector<std::vector<double>>
OrthogonalizeMatrixSpan(const std::vector<std::vector<double>>& matrix,
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

  // Create working copy of the matrix - we will orthogonalize its columns in place
  std::vector<std::vector<double>> Q = matrix;

  // Orthogonalize columns
  for (size_t j = 0; j < n; ++j)
  {
    // Two passes of orthogonalization for numerical stability
    for (int pass = 0; pass < 2; ++pass)
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
