// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include "framework/data_types/ndarray.h"
#include "framework/math/dynamic_vector.h"
#include "framework/math/dense_vector.h"

namespace opensn
{

template <typename TYPE>
class DenseMatrix : public NDArray<TYPE>
{
public:
  /// Create an empty dense matrix
  DenseMatrix() : NDArray<TYPE>() {}

  /**
   * Create a dense matrix with specified number of rows and columns
   *
   * \param rows Number of rows
   * \param cols Number of columns
   */
  DenseMatrix(unsigned int rows, unsigned int cols) : NDArray<TYPE>({rows, cols}) {}

  /**
   * Create a dense matrix with specified number of rows and columns and initialize the element to a
   * given value
   *
   * \param rows Number of rows
   * \param cols Number of columns
   * \param intitial_value Value to initialize the matrix elements with
   */
  DenseMatrix(unsigned int rows, unsigned int cols, TYPE intitial_value)
    : NDArray<TYPE>({rows, cols}, intitial_value)
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

  void SetRow(int row, const DynamicVector<TYPE>& values)
  {
    assert(Columns() == values.size());
    for (unsigned int i = 0; i < Columns(); ++i)
      (*this)(row, i) = values[i];
  }

  void SetRow(int row, const DenseVector<TYPE>& values)
  {
    assert(Columns() == values.Rows());
    for (unsigned int i = 0; i < Columns(); ++i)
      (*this)(row, i) = values(i);
  }

  /// Scale the matrix with a constant value
  void Scale(TYPE alpha)
  {
    for (unsigned int i = 0; i < Rows(); ++i)
      for (unsigned int j = 0; j < Columns(); ++j)
        (*this)(i, j) *= alpha;
  }

  /// Add matrix to this matrix
  DenseMatrix<TYPE> operator+(const DenseMatrix<TYPE>& rhs) const
  {
    assert(Rows() == rhs.Rows());
    assert(Columns() == rhs.Columns());
    DenseMatrix<TYPE> res(Rows(), Columns());
    for (unsigned int i = 0; i < Rows(); ++i)
      for (unsigned int j = 0; j < Columns(); ++j)
        res(i, j) = (*this)(i, j) + rhs(i, j);
    return res;
  }

  /// Subtract matrix from this matrix
  DenseMatrix<TYPE> operator-(const DenseMatrix<TYPE>& rhs) const
  {
    assert(Rows() == rhs.Rows());
    assert(Columns() == rhs.Columns());
    DenseMatrix<TYPE> res(Rows(), Columns());
    for (unsigned int i = 0; i < Rows(); ++i)
      for (unsigned int j = 0; j < Columns(); ++j)
        res(i, j) = (*this)(i, j) - rhs(i, j);
    return res;
  }

  /// Matrix-Vector multiplication
  DynamicVector<TYPE> operator*(const DynamicVector<TYPE>& V)
  {
    if (Rows() != V.size())
      throw std::length_error("Mismatched matrix/vector sizes in matrix-vector multiplication");

    DynamicVector<TYPE> res(Rows());
    unsigned int k = 0;
    for (unsigned int i = 0; i < Rows(); ++i)
    {
      TYPE value = 0.0;
      for (unsigned int j = 0; j < Columns(); ++j)
        value += (*this)(i, j) * V[j];
      res[k] = value;
      ++k;
    }

    return res;
  }

  /// Matrix-Vector multiplication
  DenseVector<TYPE> operator*(const DenseVector<TYPE>& b) const
  {
    if (Rows() != b.size())
      throw std::length_error("Mismatched matrix/vector sizes in matrix-vector multiplication");

    DenseVector<TYPE> res(Rows());
    unsigned int k = 0;
    for (unsigned int i = 0; i < Rows(); ++i)
    {
      TYPE value = 0.0;
      for (unsigned int j = 0; j < Columns(); ++j)
        value += (*this)(i, j) * b(j);
      res(k) = value;
      ++k;
    }

    return res;
  }

  /// Matrix-Matrix multiplication
  DenseMatrix<TYPE> operator*(const DenseMatrix<TYPE>& B) const
  {
    size_t AR = Rows();

    assert(AR != 0 and B.size() != 0);

    size_t AC = Columns();
    size_t BC = B.Columns();

    assert(AC != 0 and BC != 0 and AC == B.size());

    size_t CR = AR;
    size_t CC = BC;
    size_t Cs = AC;
    DenseMatrix<TYPE> C(CR, CC, 0.);
    for (size_t i = 0; i < CR; ++i)
      for (size_t j = 0; j < CC; ++j)
        for (size_t k = 0; k < Cs; ++k)
          C(i, j) += (*this)(i, k) * B(k, j);
    return C;
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

/// Multiplication by a scalar from the left.
template <typename TYPE>
DenseMatrix<TYPE>
operator*(double value, DenseMatrix<TYPE>& that)
{
  auto rows = that.Rows();
  auto cols = that.Columns();
  DenseMatrix<TYPE> res(rows, cols);
  for (int i = 0; i < rows; ++i)
    for (int j = 0; j < cols; ++j)
      res(i, j) = that(i, j) * value;

  return res;
}

} // namespace opensn
