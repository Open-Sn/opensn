// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include "framework/data_types/ndarray.h"
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

  void SetRow(int row, const DenseVector<TYPE>& values)
  {
    assert(Columns() == values.Rows());
    for (unsigned int i = 0; i < Columns(); ++i)
      (*this)(row, i) = values(i);
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

} // namespace opensn
