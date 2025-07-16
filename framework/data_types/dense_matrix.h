// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include "framework/data_types/ndarray.h"
#include "framework/data_types/vector.h"

namespace opensn
{

template <typename TYPE>
class DenseMatrix : public NDArray<TYPE, 2>
{
public:
  /// Create an empty dense matrix
  DenseMatrix() : NDArray<TYPE, 2>() {}

  /**
   * Create a dense matrix with specified number of rows and columns
   *
   * \param rows Number of rows
   * \param cols Number of columns
   */
  DenseMatrix(unsigned int rows, unsigned int cols) : NDArray<TYPE, 2>({rows, cols}) {}

  /**
   * Create a dense matrix with specified number of rows and columns and initialize the element to a
   * given value
   *
   * \param rows Number of rows
   * \param cols Number of columns
   * \param intitial_value Value to initialize the matrix elements with
   */
  DenseMatrix(unsigned int rows, unsigned int cols, TYPE intitial_value)
    : NDArray<TYPE, 2>({rows, cols}, intitial_value)
  {
  }

  /// Return the number of rows
  unsigned int Rows() const { return this->dimension()[0]; }

  /// Return the number of columns
  unsigned int Columns() const { return this->dimension()[1]; }

  /// Set the elements of the matrix to a specified value
  void Set(TYPE val) { this->set(val); }

  /// Set the diagonal of the matrix
  void SetDiagonal(TYPE val)
  {
    auto dims = this->dimension();
    auto d = std::min(dims[0], dims[1]);
    for (int i = 0; i < d; ++i)
      (*this)(i, i) = val;
  }

  void SetRow(int row, const Vector<TYPE>& values)
  {
    assert(Columns() == values.Rows());
    for (unsigned int i = 0; i < Columns(); ++i)
      (*this)(row, i) = values(i);
  }

  void Add(const DenseMatrix<TYPE>& other)
  {
    assert(Rows() == other.Rows());
    assert(Columns() == other.Columns());
    for (auto i = 0; i < Rows(); ++i)
      for (auto j = 0; j < Columns(); ++j)
        (*this)(i, j) += other(i, j);
  }

  void Subtract(const DenseMatrix<TYPE>& other)
  {
    assert(Rows() == other.Rows());
    assert(Columns() == other.Columns());
    for (auto i = 0; i < Rows(); ++i)
      for (auto j = 0; j < Columns(); ++j)
        (*this)(i, j) -= other(i, j);
  }

  /// Multiply matrix with a vector and return resulting vector
  Vector<TYPE> Mult(const Vector<TYPE>& x) const
  {
    auto rows = Rows();
    auto cols = x.Rows();

    assert(rows > 0);
    assert(cols == Columns());

    Vector<TYPE> b(rows, 0.0);
    for (auto i = 0; i < rows; ++i)
      for (auto j = 0; j < cols; ++j)
        b(i) += (*this)(i, j) * x(j);
    return b;
  }

  /// Multiply by a matrix and return the resulting matrix
  DenseMatrix<TYPE> Mult(const DenseMatrix<TYPE>& other) const
  {
    auto rows = Rows();
    assert(rows != 0 and other.Rows() != 0);
    auto columns = Columns();
    auto other_cols = other.Columns();
    assert(columns != 0 and other_cols != 0 and columns == other.Rows());
    DenseMatrix<TYPE> res(rows, other_cols, 0.);
    for (auto i = 0; i < rows; ++i)
      for (auto j = 0; j < other_cols; ++j)
        for (auto k = 0; k < columns; ++k)
          res(i, j) += (*this)(i, k) * other(k, j);
    return res;
  }

  /// Returns the transpose of a matrix
  DenseMatrix<TYPE> Transposed() const
  {
    assert(Rows());
    assert(Columns());
    DenseMatrix<TYPE> trans(Columns(), Rows());
    for (auto i = 0; i < Rows(); ++i)
      for (auto j = 0; j < Columns(); ++j)
        trans(j, i) = (*this)(i, j);
    return trans;
  }

  /// Transpose this matrix
  void Transpose()
  {
    assert(Rows());
    assert(Columns());
    DenseMatrix<TYPE> trans(Columns(), Rows());
    for (auto i = 0; i < Rows(); ++i)
      for (auto j = 0; j < Columns(); ++j)
        trans(j, i) = (*this)(i, j);
    (*this) = trans;
  }

  /// Prints the matrix to a string and then returns the string.
  std::string PrintStr() const
  {
    std::stringstream out;

    for (int i = 0; i < Rows(); ++i)
    {
      for (int j = 0; j < (Columns() - 1); ++j)
        out << (*this)(i, j) << " ";
      out << (*this)(i, Columns() - 1);

      if (i < (Rows() - 1))
        out << "\n";
    }

    return out.str();
  }
};

/// Returns the transpose of a matrix.
template <typename TYPE>
DenseMatrix<TYPE>
Transpose(const DenseMatrix<TYPE>& A)
{
  assert(A.Rows());
  assert(A.Columns());
  auto AR = A.Rows();
  auto AC = 0;
  if (AR)
    AC = A.Columns();

  DenseMatrix<TYPE> T(AC, AR);
  for (auto i = 0; i < AR; ++i)
    for (auto j = 0; j < AC; ++j)
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

  for (auto j = 0; j < cols; ++j)
    std::swap(A(r1, j), A(r2, j));
}

/// Multiply matrix with a constant and return result.
template <typename TYPE>
DenseMatrix<TYPE>
Mult(const DenseMatrix<TYPE>& A, const TYPE c)
{
  auto R = A.Rows();
  auto C = A.Columns();
  DenseMatrix<TYPE> B(R, C, 0.);
  for (auto i = 0; i < R; ++i)
    for (auto j = 0; j < C; ++j)
      B(i, j) = A(i, j) * c;
  return B;
}

/// Multiply matrix with a vector and return resulting vector
template <typename TYPE>
Vector<TYPE>
Mult(const DenseMatrix<TYPE>& A, const Vector<TYPE>& x)
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
Mult(const DenseMatrix<TYPE>& A, const DenseMatrix<TYPE>& B)
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
Add(const DenseMatrix<TYPE>& A, const DenseMatrix<TYPE>& B)
{
  auto AR = A.Rows();
  auto BR = B.Rows();

  assert(AR != 0 and B.Rows() != 0);
  assert(AR == BR);

  auto AC = A.Columns();
  auto BC = B.Columns();

  assert(AC != 0 and BC != 0);
  assert(AC == BC);

  DenseMatrix<TYPE> C(AR, AC, 0.0);
  for (auto i = 0; i < AR; ++i)
    for (auto j = 0; j < AC; ++j)
      C(i, j) = A(i, j) + B(i, j);
  return C;
}

/// Subtracts matrix A from B and returns the result.
template <typename TYPE>
DenseMatrix<TYPE>
Subtract(const DenseMatrix<TYPE>& A, const DenseMatrix<TYPE>& B)
{
  auto AR = A.Rows();
  auto BR = B.Rows();

  assert(AR != 0 and B.size() != 0);
  assert(AR == BR);

  auto AC = A.Columns();
  auto BC = B.Columns();

  assert(AC != 0 and BC != 0);
  assert(AC == BC);

  DenseMatrix<TYPE> C(AR, AC, 0.0);
  for (auto i = 0; i < AR; ++i)
    for (auto j = 0; j < AC; ++j)
      C(i, j) = A(i, j) - B(i, j);
  return C;
}

/// Scale the matrix with a constant value
template <typename TYPE>
void
Scale(DenseMatrix<TYPE>& mat, TYPE alpha)
{
  for (auto i = 0; i < mat.Rows(); ++i)
    for (auto j = 0; j < mat.Columns(); ++j)
      mat(i, j) *= alpha;
}

/// Scale the matrix with a constant value
template <typename TYPE>
DenseMatrix<TYPE>
Scaled(const DenseMatrix<TYPE>& mat, TYPE alpha)
{
  DenseMatrix<TYPE> res(mat.Rows(), mat.Columns());
  for (auto i = 0; i < mat.Rows(); ++i)
    for (auto j = 0; j < mat.Columns(); ++j)
      res(i, j) = mat(i, j) * alpha;
  return res;
}

/// Returns a copy of A with removed rows `r` and removed columns `c`
template <typename TYPE>
DenseMatrix<TYPE>
SubMatrix(const DenseMatrix<TYPE>& A, const size_t r, const size_t c)
{
  auto rows = A.Rows();
  auto cols = A.Columns();
  assert((r >= 0) and (r < rows) and (c >= 0) and (c < cols));

  DenseMatrix<TYPE> B(rows - 1, cols - 1);
  for (auto i = 0, ii = 0; i < rows; ++i)
  {
    if (i != r)
    {
      for (auto j = 0, jj = 0; j < cols; ++j)
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
  auto rows = A.Rows();

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
    for (auto n = 0; n < rows; ++n)
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
  for (auto i = 0; i < n - 1; ++i)
  {
    auto bi = b(i);
    auto factor = 1.0 / A(i, i);
    for (auto j = i + 1; j < n; ++j)
    {
      auto val = A(j, i) * factor;
      b(j) -= val * bi;
      for (auto k = i + 1; k < n; ++k)
        A(j, k) -= val * A(i, k);
    }
  }

  // Back substitution
  for (int i = n - 1; i >= 0; --i)
  {
    auto bi = b(i);
    for (auto j = i + 1; j < n; ++j)
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

  const auto rows = A.Rows();
  const auto cols = A.Columns();

  DenseMatrix<TYPE> M(rows, cols);
  M.Set(0.);
  M.SetDiagonal(1.);

  auto B = A;

  for (auto c = 0; c < rows; ++c)
  {
    // Find a row with the largest pivot value
    auto max_row = c; // nzr = non-zero row
    for (auto r = c; r < rows; ++r)
      if (std::fabs(B(r, c)) > std::fabs(B(max_row, c)))
        max_row = r;

    if (max_row != c)
    {
      SwapRows(B, max_row, c);
      SwapRows(M, max_row, c);
    }

    // Eliminate non-zero values
    for (auto r = 0; r < rows; ++r)
    {
      if (r != c)
      {
        auto g = B(r, c) / B(c, c);
        if (B(r, c) != 0)
        {
          for (auto k = 0; k < rows; ++k)
          {
            B(r, k) -= B(c, k) * g;
            M(r, k) -= M(c, k) * g;
          }
          B(r, c) = 0;
        }
      }
      else
      {
        auto g = 1 / B(c, c);
        for (auto k = 0; k < rows; ++k)
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
  auto rows = A.Rows();
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
  auto n = A.Rows();
  int it_counter = 0;
  Vector<double> y(n, 1.0);
  double lambda0 = 0.0;

  // Perform initial iteration outside of loop
  auto Ay = Mult(A, y);
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
    Ay = Mult(A, y);
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

} // namespace opensn
