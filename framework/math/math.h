// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include "framework/math/math_incdef.h"
#include "framework/math/quadratures/spatial/spatial_quadrature.h"
#include "framework/math/quadratures/angular/angular_quadrature.h"
#include "framework/math/unknown_manager/unknown_manager.h"
#include "framework/math/dense_vector.h"
#include "framework/math/dense_matrix.h"
#include <memory>

namespace opensn
{
class SparseMatrix;
class UnknownManager;
class CDFSampler;
class SpatialDiscretization;
class SpatialDiscretization_FV;
class SpatialDiscretization_PWLD;
class SpatialDiscretization_PWLC;

using MatVec3 = std::vector<std::vector<Vector3>>;

/// Coordinate system type.
enum class CoordinateSystemType
{
  UNDEFINED = 0,
  CARTESIAN = 1,
  CYLINDRICAL = 2,
  SPHERICAL = 3,
};

/// Spatial discretization type.
enum class SpatialDiscretizationType
{
  UNDEFINED = 0,
  FINITE_VOLUME = 1,
  PIECEWISE_LINEAR_CONTINUOUS = 2,
  PIECEWISE_LINEAR_DISCONTINUOUS = 3,
  LAGRANGE_CONTINUOUS = 4,
  LAGRANGE_DISCONTINUOUS = 5
};

enum class NormType : int
{
  L1_NORM = 1,
  L2_NORM = 2,
  LINF_NORM = 3
};

/**
 * Sample a Cumulative Distribution Function (CDF) given a probability.
 *
 * The supplied vector should contain the upper bin boundary for each
 * bin and will return the bin associated with the bin that brackets
 * the supplied probability.
 *
 * Example:
 * Suppose we sample bins 0-9. Suppose also that the probalities for each
 * bin is as follows:
 * - 0.1 bin 0
 * - 0.1 bin 1
 * - 0.5 bin 5
 * - 0.3 bin 8
 *
 * The CDF for this probability distribution will look like this
 * - bin 0 = 0.1
 * - bin 1 = 0.2
 * - bin 2 = 0.2
 * - bin 3 = 0.2
 * - bin 4 = 0.2
 * - bin 5 = 0.7
 * - bin 6 = 0.7
 * - bin 7 = 0.7
 * - bin 8 = 1.0
 * - bin 9 = 1.0
 *
 * Supplying a random number between 0 and 1 should indicate sampling one
 * of the bins 0,1,5 or 8. The most inefficient way to do this is to
 * linearly loop through the cdf and check \f$ cdf_{i-1} \ge \theta < cdf_i \f$.
 *  An optimized version of this sampling would be to perform a recursive
 *  block search which starts with a course view of the cdf and then gradually
 *  refines the view until the final linear search can be performed.*/
int SampleCDF(double x, std::vector<double> cdf_bin);

/// Computes the factorial of an integer.
double Factorial(int x);

/**
 * Determines the azimuthal- and polar-angle associated with the given direction vector.
 * Returns a pair = [azimuthal-angle,polar-angle].
 */
std::pair<double, double> OmegaToPhiThetaSafe(const Vector3& omega);

/// Prints the Vector.
void PrintVector(const std::vector<double>& x);

/// Scales a vector in place by constant.
void Scale(std::vector<double>& x, const double& val);

/// Scale the vector with a constant value
template <typename TYPE>
void
Scale(Vector<TYPE>& a, TYPE alpha)
{
  for (unsigned int i = 0; i < a.Rows(); ++i)
    a(i) *= alpha;
}

template <typename TYPE>
Vector<TYPE>
Scaled(const Vector<TYPE>& a, TYPE alpha)
{
  Vector<TYPE> res(a.Rows());
  for (unsigned int i = 0; i < a.Rows(); ++i)
    res(i) = a(i) * alpha;
  return res;
}

/// Sets a constant value to a vector.
void Set(std::vector<double>& x, const double& val);

/// Add vector to this vector
template <typename TYPE>
Vector<TYPE>
VecAdd(const Vector<TYPE>& a, const Vector<TYPE>& b)
{
  assert(a.Rows() == b.Rows());
  Vector<TYPE> res(a.Rows());
  for (unsigned int i = 0; i < a.Rows(); ++i)
    res(i) = a(i) + b(i);
  return res;
}

/// Subtract two vectors
template <typename TYPE>
Vector<TYPE>
VecSub(const Vector<TYPE>& a, const Vector<TYPE>& b)
{
  assert(a.Rows() == b.Rows());
  Vector<TYPE> res(a.Rows());
  for (unsigned int i = 0; i < a.Rows(); ++i)
    res(i) = a(i) - b(i);
  return res;
}

/// Multiplies the vector with a constant and returns result.
std::vector<double> VecMul(const std::vector<double>& x, const double& val);

/**
 * Returns the 1-norm. Also known as the Taxicab or Manhattan norm.
 *
 * \f[
 * \|\boldsymbol{x}\|_{1}=\sum_{i=1}^{n}\left|x_{i}\right|
 * \f]
 */
double Vec1Norm(const std::vector<double>& x);

/**
 * Returns the 2-norm. Also known as the Euclidian or Frobenius norm.
 *
 * \f[
 * \|\boldsymbol{x}\|_{2}=\sqrt{x_{1}^{2}+\cdots+x_{n}^{2}}
 * \f]
 */
double Vec2Norm(const std::vector<double>& x);

template <typename TYPE>
double
Vec2Norm(const Vector<TYPE>& x)
{
  size_t n = x.Rows();
  double val = 0.0;
  for (size_t i = 0; i != n; ++i)
    val += x(i) * x(i);
  return std::sqrt(val);
}

/**
 * Returns the infinity-norm.
 *
 * \f[
 * \|\mathbf{x}\|_{\infty}=\max \left(\left|x_{1}\right|,
 * \ldots,\left|x_{n}\right|\right) \f]
 */
double VecInfinityNorm(const std::vector<double>& x);

/**
 * Returns the p-norm.
 *
 * \f[
 * \|\mathbf{x}\|_{p}=\left(\sum_{i=1}^{n}\left|x_{i}\right|^{p}\right)^{1 / p}
 * \f]
 */
double VecPNorm(const std::vector<double>& x, const double& p);

/**
 * Computes the dot product of two vectors.
 *
 * \f[
 * \mathrm{a} \cdot \mathrm{b}=\sum_{i=1}^{n} a_{i} b_{i}
 * \f]
 */
double Dot(const std::vector<double>& x, const std::vector<double>& y);

template <typename TYPE>
double
Dot(const Vector<TYPE>& x, const Vector<TYPE>& y)
{
  assert(x.Rows() > 0);
  assert(y.Rows() > 0);
  assert(x.Rows() == y.Rows());
  double val = 0.0;
  for (size_t i = 0; i < x.Rows(); ++i)
    val += x(i) * y(i);
  return val;
}

/// Adds two vectors component-wise.
std::vector<double> operator+(const std::vector<double>& a, const std::vector<double>& b);

/// Subtracts two vectors component-wise.
std::vector<double> operator-(const std::vector<double>& a, const std::vector<double>& b);

/// Returns the transpose of a matrix.
template <typename TYPE>
DenseMatrix<TYPE>
Transpose(const DenseMatrix<TYPE>& A)
{
  assert(A.Rows());
  assert(A.Column());
  size_t AR = A.Rows();
  size_t AC = 0;
  if (AR)
    AC = A.Columns();

  DenseMatrix<TYPE> T(AC, AR);
  for (size_t i = 0; i < AR; ++i)
    for (size_t j = 0; j < AC; ++j)
      T(j, i) = A(i, j);
  return T;
}

/// Swaps two rows of a matrix.
template <typename TYPE>
void
SwapRows(DenseMatrix<TYPE>& A, size_t r1, size_t r2)
{
  auto rows = A.Rows();
  auto cols = A.Columns();
  assert(rows > 0);
  assert(cols > 0);
  assert(r1 >= 0 and r1 < rows and r2 >= 0 and r2 < rows);

  for (size_t j = 0; j < cols; ++j)
    std::swap(A(r1, j), A(r2, j));
}

/// Multiply matrix with a constant and return result.
template <typename TYPE>
DenseMatrix<TYPE>
MatMul(const DenseMatrix<TYPE>& A, const TYPE c)
{
  auto R = A.Rows();
  auto C = A.Rows();
  DenseMatrix<TYPE> B(R, C, 0.);
  for (auto i = 0; i < R; ++i)
    for (auto j = 0; j < C; ++j)
      B(i, j) = A(i, j) * c;
  return B;
}

/// Multiply matrix with a vector and return resulting vector
template <typename TYPE>
Vector<TYPE>
MatMul(const DenseMatrix<TYPE>& A, const Vector<TYPE>& x)
{
  auto R = A.Rows();
  auto C = x.Rows();

  assert(R > 0);
  assert(C == A.Columns());

  Vector<TYPE> b(R, 0.0);
  for (auto i = 0; i < R; ++i)
  {
    for (auto j = 0; j < C; ++j)
      b(i) += A(i, j) * x(j);
  }

  return b;
}

/// Mutliply two matrices and return result.
template <typename TYPE>
DenseMatrix<TYPE>
MatMul(const DenseMatrix<TYPE>& A, const DenseMatrix<TYPE>& B)
{
  auto AR = A.Rows();

  assert(AR != 0 and B.Rows() != 0);

  auto AC = A.Columns();
  auto BC = B.Columns();

  assert(AC != 0 and BC != 0 and AC == B.Rows());

  auto CR = AR;
  auto CC = BC;
  auto Cs = AC;

  DenseMatrix<TYPE> C(CR, CC, 0.);
  for (auto i = 0; i < CR; ++i)
    for (auto j = 0; j < CC; ++j)
      for (auto k = 0; k < Cs; ++k)
        C(i, j) += A(i, k) * B(k, j);
  return C;
}

/// Adds two matrices and returns the result.
template <typename TYPE>
DenseMatrix<TYPE>
MatAdd(const DenseMatrix<TYPE>& A, const DenseMatrix<TYPE>& B)
{
  size_t AR = A.Rows();
  size_t BR = B.Rows();

  assert(AR != 0 and B.Rows() != 0);
  assert(AR == BR);

  size_t AC = A.Columns();
  size_t BC = B.Columns();

  assert(AC != 0 and BC != 0);
  assert(AC == BC);

  DenseMatrix<TYPE> C(AR, AC, 0.0);
  for (size_t i = 0; i < AR; ++i)
    for (size_t j = 0; j < AC; ++j)
      C(i, j) = A(i, j) + B(i, j);
  return C;
}

/// Subtracts matrix A from B and returns the result.
template <typename TYPE>
DenseMatrix<TYPE>
MatSubtract(const DenseMatrix<TYPE>& A, const DenseMatrix<TYPE>& B)
{
  size_t AR = A.Rows();
  size_t BR = B.Rows();

  assert(AR != 0 and B.size() != 0);
  assert(AR == BR);

  size_t AC = A.Columns();
  size_t BC = B.Columns();

  assert(AC != 0 and BC != 0);
  assert(AC == BC);

  DenseMatrix<TYPE> C(AR, AC, 0.0);
  for (size_t i = 0; i < AR; ++i)
    for (size_t j = 0; j < AC; ++j)
      C(i, j) = A(i, j) - B(i, j);
  return C;
}

/// Scale the matrix with a constant value
template <typename TYPE>
void
Scale(DenseMatrix<TYPE>& mat, TYPE alpha)
{
  for (unsigned int i = 0; i < mat.Rows(); ++i)
    for (unsigned int j = 0; j < mat.Columns(); ++j)
      mat(i, j) *= alpha;
}

/// Scale the matrix with a constant value
template <typename TYPE>
DenseMatrix<TYPE>
Scaled(const DenseMatrix<TYPE>& mat, TYPE alpha)
{
  DenseMatrix<TYPE> res(mat.Rows(), mat.Columns());
  for (unsigned int i = 0; i < mat.Rows(); ++i)
    for (unsigned int j = 0; j < mat.Columns(); ++j)
      res(i, j) = mat(i, j) * alpha;
  return res;
}

/// Returns a copy of A with removed rows `r` and removed columns `c`
template <typename TYPE>
DenseMatrix<TYPE>
SubMatrix(const DenseMatrix<TYPE>& A, const size_t r, const size_t c)
{
  size_t rows = A.Rows();
  size_t cols = A.Columns();
  assert((r >= 0) and (r < rows) and (c >= 0) and (c < cols));

  DenseMatrix<TYPE> B(rows - 1, cols - 1);
  for (size_t i = 0, ii = 0; i < rows; ++i)
  {
    if (i != r)
    {
      for (size_t j = 0, jj = 0; j < cols; ++j)
      {
        if (j != c)
        {
          B(ii, jj) = A(i, j);
          ++jj;
        }
      }
      ++ii;
    }
  }
  return B;
}

/// Computes the determinant of a matrix.
template <typename TYPE>
double
Determinant(const DenseMatrix<TYPE>& A)
{
  size_t rows = A.Rows();

  if (rows == 1)
    return A(0, 0);
  else if (rows == 2)
  {
    return A(0, 0) * A(1, 1) - A(0, 1) * A(1, 0);
  }
  else if (rows == 3)
  {
    return A(0, 0) * A(1, 1) * A(2, 2) + A(0, 1) * A(1, 2) * A(2, 0) + A(0, 2) * A(1, 0) * A(2, 1) -
           A(0, 0) * A(1, 2) * A(2, 1) - A(0, 1) * A(1, 0) * A(2, 2) - A(0, 2) * A(1, 1) * A(2, 0);
  }
  // http://www.cvl.iis.u-tokyo.ac.jp/~Aiyazaki/tech/teche23.htAl
  else if (rows == 4)
  {
    return A(0, 0) * A(1, 1) * A(2, 2) * A(3, 3) + A(0, 0) * A(1, 2) * A(2, 3) * A(3, 1) +
           A(0, 0) * A(1, 3) * A(2, 1) * A(3, 2) + A(0, 1) * A(1, 0) * A(2, 3) * A(3, 2) +
           A(0, 1) * A(1, 2) * A(2, 0) * A(3, 3) + A(0, 1) * A(1, 3) * A(2, 2) * A(3, 0) +
           A(0, 2) * A(1, 0) * A(2, 1) * A(3, 3) + A(0, 2) * A(1, 1) * A(2, 3) * A(3, 0) +
           A(0, 2) * A(1, 3) * A(2, 0) * A(3, 1) + A(0, 3) * A(1, 0) * A(2, 2) * A(3, 1) +
           A(0, 3) * A(1, 1) * A(2, 0) * A(3, 2) + A(0, 3) * A(1, 2) * A(2, 1) * A(3, 0) -
           A(0, 0) * A(1, 1) * A(2, 3) * A(3, 2) - A(0, 0) * A(1, 2) * A(2, 1) * A(3, 3) -
           A(0, 0) * A(1, 3) * A(2, 2) * A(3, 1) - A(0, 1) * A(1, 0) * A(2, 2) * A(3, 3) -
           A(0, 1) * A(1, 2) * A(2, 3) * A(3, 0) - A(0, 1) * A(1, 3) * A(2, 0) * A(3, 2) -
           A(0, 2) * A(1, 0) * A(2, 3) * A(3, 1) - A(0, 2) * A(1, 1) * A(2, 0) * A(3, 3) -
           A(0, 2) * A(1, 3) * A(2, 1) * A(3, 0) - A(0, 3) * A(1, 0) * A(2, 1) * A(3, 2) -
           A(0, 3) * A(1, 1) * A(2, 2) * A(3, 0) - A(0, 3) * A(1, 2) * A(2, 0) * A(3, 1);
  }
  else
  {
    double det = 0;
    for (size_t n = 0; n < rows; ++n)
    {
      auto M = SubMatrix(A, 0, n);
      double pm = ((n + 1) % 2) * 2.0 - 1.0;
      det += pm * A(0, n) * Determinant(M);
    }
    return det;
  }
}

/// Gauss Elimination without pivoting.
template <typename TYPE>
void
GaussElimination(DenseMatrix<TYPE>& A, Vector<TYPE>& b, unsigned int n)
{
  // Forward elimination
  for (unsigned int i = 0; i < n - 1; ++i)
  {
    auto bi = b(i);
    auto factor = 1.0 / A(i, i);
    for (unsigned int j = i + 1; j < n; ++j)
    {
      auto val = A(j, i) * factor;
      b(j) -= val * bi;
      for (unsigned int k = i + 1; k < n; ++k)
        A(j, k) -= val * A(i, k);
    }
  }

  // Back substitution
  for (int i = n - 1; i >= 0; --i)
  {
    auto bi = b(i);
    for (unsigned int j = i + 1; j < n; ++j)
      bi -= A(i, j) * b(j);
    b(i) = bi / A(i, i);
  }
}

/// Computes the inverse of a matrix using Gauss-Elimination with pivoting.
template <typename TYPE>
DenseMatrix<TYPE>
InverseGEPivoting(const DenseMatrix<TYPE>& A)
{
  assert(A.Rows() == A.Columns());

  const unsigned int rows = A.Rows();
  const unsigned int cols = A.Columns();

  DenseMatrix<TYPE> M(rows, cols);
  M.Set(0.);
  M.SetDiagonal(1.);

  auto B = A;

  for (unsigned int c = 0; c < rows; ++c)
  {
    // Find a row with the largest pivot value
    unsigned int max_row = c; // nzr = non-zero row
    for (unsigned int r = c; r < rows; ++r)
      if (std::fabs(B(r, c)) > std::fabs(B(max_row, c)))
        max_row = r;

    if (max_row != c)
    {
      SwapRows(B, max_row, c);
      SwapRows(M, max_row, c);
    }

    // Eliminate non-zero values
    for (unsigned int r = 0; r < rows; ++r)
    {
      if (r != c)
      {
        double g = B(r, c) / B(c, c);
        if (B(r, c) != 0)
        {
          for (unsigned int k = 0; k < rows; ++k)
          {
            B(r, k) -= B(c, k) * g;
            M(r, k) -= M(c, k) * g;
          }
          B(r, c) = 0;
        }
      }
      else
      {
        double g = 1 / B(c, c);
        for (unsigned int k = 0; k < rows; ++k)
        {
          B(r, k) *= g;
          M(r, k) *= g;
        }
      }
    }
  }
  return M;
}

/// Computes the inverse of a matrix.
template <typename TYPE>
DenseMatrix<TYPE>
Inverse(const DenseMatrix<TYPE>& A)
{
  size_t rows = A.Rows();
  DenseMatrix<double> M(rows, A.Rows());
  double f = 0.0;

  // Only calculate the determinant if matrix size is less than 5 since
  // the inverse is directly calculated for larger matrices. Otherwise,
  // the inverse routine spends all of its time sitting in the determinant
  // function which is unnecessary.
  if (rows < 5)
  {
    f = Determinant(A);
    assert(f != 0.0);
    f = 1.0 / f;
  }

  if (rows == 1)
    M(0, 0) = f;
  else if (rows == 2)
  {
    M(0, 0) = A(1, 1);
    M(0, 1) = -A(0, 1);
    M(1, 0) = -A(1, 0);
    M(1, 1) = A(0, 0);
    Scale(M, f);
  }
  else if (rows == 3)
  {
    M(0, 0) = A(2, 2) * A(1, 1) - A(2, 1) * A(1, 2);
    M(0, 1) = -(A(2, 2) * A(0, 1) - A(2, 1) * A(0, 2));
    M(0, 2) = A(1, 2) * A(0, 1) - A(1, 1) * A(0, 2);
    M(1, 0) = -(A(2, 2) * A(1, 0) - A(2, 0) * A(1, 2));
    M(1, 1) = A(2, 2) * A(0, 0) - A(2, 0) * A(0, 2);
    M(1, 2) = -(A(1, 2) * A(0, 0) - A(1, 0) * A(0, 2));
    M(2, 0) = A(2, 1) * A(1, 0) - A(2, 0) * A(1, 1);
    M(2, 1) = -(A(2, 1) * A(0, 0) - A(2, 0) * A(0, 1));
    M(2, 2) = A(1, 1) * A(0, 0) - A(1, 0) * A(0, 1);
    Scale(M, f);
  }
  else if (rows == 4)
  {
    // http://www.cvl.iis.u-tokyo.ac.jp/~Aiyazaki/tech/teche23.htAl
    M(0, 0) = A(1, 1) * A(2, 2) * A(3, 3) + A(1, 2) * A(2, 3) * A(3, 1) +
              A(1, 3) * A(2, 1) * A(3, 2) - A(1, 1) * A(2, 3) * A(3, 2) -
              A(1, 2) * A(2, 1) * A(3, 3) - A(1, 3) * A(2, 2) * A(3, 1);
    M(0, 1) = A(0, 1) * A(2, 3) * A(3, 2) + A(0, 2) * A(2, 1) * A(3, 3) +
              A(0, 3) * A(2, 2) * A(3, 1) - A(0, 1) * A(2, 2) * A(3, 3) -
              A(0, 2) * A(2, 3) * A(3, 1) - A(0, 3) * A(2, 1) * A(3, 2);
    M(0, 2) = A(0, 1) * A(1, 2) * A(3, 3) + A(0, 2) * A(1, 3) * A(3, 1) +
              A(0, 3) * A(1, 1) * A(3, 2) - A(0, 1) * A(1, 3) * A(3, 2) -
              A(0, 2) * A(1, 1) * A(3, 3) - A(0, 3) * A(1, 2) * A(3, 1);
    M(0, 3) = A(0, 1) * A(1, 3) * A(2, 2) + A(0, 2) * A(1, 1) * A(2, 3) +
              A(0, 3) * A(1, 2) * A(2, 1) - A(0, 1) * A(1, 2) * A(2, 3) -
              A(0, 2) * A(1, 3) * A(2, 1) - A(0, 3) * A(1, 1) * A(2, 2);

    M(1, 0) = A(1, 0) * A(2, 3) * A(3, 2) + A(1, 2) * A(2, 0) * A(3, 3) +
              A(1, 3) * A(2, 2) * A(3, 0) - A(1, 0) * A(2, 2) * A(3, 3) -
              A(1, 2) * A(2, 3) * A(3, 0) - A(1, 3) * A(2, 0) * A(3, 2);
    M(1, 1) = A(0, 0) * A(2, 2) * A(3, 3) + A(0, 2) * A(2, 3) * A(3, 0) +
              A(0, 3) * A(2, 0) * A(3, 2) - A(0, 0) * A(2, 3) * A(3, 2) -
              A(0, 2) * A(2, 0) * A(3, 3) - A(0, 3) * A(2, 2) * A(3, 0);
    M(1, 2) = A(0, 0) * A(1, 3) * A(3, 2) + A(0, 2) * A(1, 0) * A(3, 3) +
              A(0, 3) * A(1, 2) * A(3, 0) - A(0, 0) * A(1, 2) * A(3, 3) -
              A(0, 2) * A(1, 3) * A(3, 0) - A(0, 3) * A(1, 0) * A(3, 2);
    M(1, 3) = A(0, 0) * A(1, 2) * A(2, 3) + A(0, 2) * A(1, 3) * A(2, 0) +
              A(0, 3) * A(1, 0) * A(2, 2) - A(0, 0) * A(1, 3) * A(2, 2) -
              A(0, 2) * A(1, 0) * A(2, 3) - A(0, 3) * A(1, 2) * A(2, 0);

    M(2, 0) = A(1, 0) * A(2, 1) * A(3, 3) + A(1, 1) * A(2, 3) * A(3, 0) +
              A(1, 3) * A(2, 0) * A(3, 1) - A(1, 0) * A(2, 3) * A(3, 1) -
              A(1, 1) * A(2, 0) * A(3, 3) - A(1, 3) * A(2, 1) * A(3, 0);
    M(2, 1) = A(0, 0) * A(2, 3) * A(3, 1) + A(0, 1) * A(2, 0) * A(3, 3) +
              A(0, 3) * A(2, 1) * A(3, 0) - A(0, 0) * A(2, 1) * A(3, 3) -
              A(0, 1) * A(2, 3) * A(3, 0) - A(0, 3) * A(2, 0) * A(3, 1);
    M(2, 2) = A(0, 0) * A(1, 1) * A(3, 3) + A(0, 1) * A(1, 3) * A(3, 0) +
              A(0, 3) * A(1, 0) * A(3, 1) - A(0, 0) * A(1, 3) * A(3, 1) -
              A(0, 1) * A(1, 0) * A(3, 3) - A(0, 3) * A(1, 1) * A(3, 0);
    M(2, 3) = A(0, 0) * A(1, 3) * A(2, 1) + A(0, 1) * A(1, 0) * A(2, 3) +
              A(0, 3) * A(1, 1) * A(2, 0) - A(0, 0) * A(1, 1) * A(2, 3) -
              A(0, 1) * A(1, 3) * A(2, 0) - A(0, 3) * A(1, 0) * A(2, 1);

    M(3, 0) = A(1, 0) * A(2, 2) * A(3, 1) + A(1, 1) * A(2, 0) * A(3, 2) +
              A(1, 2) * A(2, 1) * A(3, 0) - A(1, 0) * A(2, 1) * A(3, 2) -
              A(1, 1) * A(2, 2) * A(3, 0) - A(1, 2) * A(2, 0) * A(3, 1);
    M(3, 1) = A(0, 0) * A(2, 1) * A(3, 2) + A(0, 1) * A(2, 2) * A(3, 0) +
              A(0, 2) * A(2, 0) * A(3, 1) - A(0, 0) * A(2, 2) * A(3, 1) -
              A(0, 1) * A(2, 0) * A(3, 2) - A(0, 2) * A(2, 1) * A(3, 0);
    M(3, 2) = A(0, 0) * A(1, 2) * A(3, 1) + A(0, 1) * A(1, 0) * A(3, 2) +
              A(0, 2) * A(1, 1) * A(3, 0) - A(0, 0) * A(1, 1) * A(3, 2) -
              A(0, 1) * A(1, 2) * A(3, 0) - A(0, 2) * A(1, 0) * A(3, 1);
    M(3, 3) = A(0, 0) * A(1, 1) * A(2, 2) + A(0, 1) * A(1, 2) * A(2, 0) +
              A(0, 2) * A(1, 0) * A(2, 1) - A(0, 0) * A(1, 2) * A(2, 1) -
              A(0, 1) * A(1, 0) * A(2, 2) - A(0, 2) * A(1, 1) * A(2, 0);
    Scale(M, f);
  }
  else
    M = InverseGEPivoting(A);

  return M;
}

/**
 * Performs power iteration to obtain the fundamental eigen mode. The eigen-value of the fundamental
 * mode is return whilst the eigen-vector is return via reference.
 */
template <typename TYPE>
double
PowerIteration(const DenseMatrix<TYPE>& A,
               Vector<TYPE>& e_vec,
               int max_it = 2000,
               double tol = 1.0e-13)
{
  unsigned int n = A.Rows();
  int it_counter = 0;
  Vector<double> y(n, 1.0);
  double lambda0 = 0.0;

  // Perform initial iteration outside of loop
  auto Ay = MatMul(A, y);
  auto lambda = Dot(y, Ay);
  y = Scaled(Ay, 1.0 / Vec2Norm(Ay));
  if (lambda < 0.0)
    Scale(y, -1.0);

  // Perform convergence loop
  bool converged = false;
  while (not converged and it_counter < max_it)
  {
    // Update old eigenvalue
    lambda0 = std::fabs(lambda);
    // Calculate new eigenvalue/eigenvector
    Ay = MatMul(A, y);
    lambda = Dot(y, Ay);
    y = Scaled(Ay, 1.0 / Vec2Norm(Ay));

    // Check if converged or not
    if (std::fabs(std::fabs(lambda) - lambda0) <= tol)
      converged = true;
    // Update counter
    ++it_counter;
  }

  if (lambda < 0.0)
    Scale(y, -1.0);

  // Renormalize eigenvector for the last time
  y = Scaled(Ay, 1.0 / lambda);

  // Set eigenvector, return the eigenvalue
  e_vec = std::move(y);

  return lambda;
}

double ComputePointwiseChange(std::vector<double>& x, std::vector<double>& y);

double ComputeL2Change(std::vector<double>& x, std::vector<double>& y);

} // namespace opensn
