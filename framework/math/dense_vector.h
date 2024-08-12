// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include "framework/data_types/ndarray.h"
#include <sstream>

namespace opensn
{

template <typename TYPE>
class DenseVector : public NDArray<TYPE>
{
public:
  /// Create an empty dense (column) vector
  DenseVector() : NDArray<TYPE>() {}

  /// Create a dense (column) vector with specified number of rows
  DenseVector(unsigned int rows) : NDArray<TYPE>({rows}) {}

  /// Create a dense (column) vector with specified number of rows and initialize the elements
  DenseVector(unsigned int rows, TYPE intitial_value) : NDArray<TYPE>({rows}, intitial_value) {}

  DenseVector(const std::vector<TYPE>& in) : NDArray<TYPE>({in.size()})
  {
    for (std::size_t i = 0; i < in.size(); ++i)
      (*this)(i) = in[i];
  }

  /// Return the number of rows
  unsigned int Rows() const
  {
    if (this->empty())
      return 0;
    else
      return this->dimension()[0];
  }

  /// Resize the vector
  void Resize(unsigned int new_size) { this->resize({new_size}); }

  /// Resize the vector and initialize the new elements with a value
  void Resize(unsigned int new_size, const TYPE& value)
  {
    auto old_size = Rows();
    this->resize({new_size});
    if (Rows() > old_size)
    {
      for (unsigned int i = old_size; i < Rows(); ++i)
        (*this)(i) = value;
    }
  }

  /// Set the elements of the vector to a specified value
  void Set(TYPE val) { this->set(val); }

  /// Scale the vector with a constant value
  void Scale(TYPE alpha)
  {
    for (unsigned int i = 0; i < Rows(); ++i)
      (*this)(i) *= alpha;
  }

  /// Normalizes the vector in-place.
  void Normalize()
  {
    TYPE norm = this->Norm();
    for (unsigned int i = 0; i < Rows(); ++i)
      (*this)(i) /= norm;
  }

  /// Add vector to this vector
  DenseVector<TYPE> operator+(const DenseVector<TYPE>& rhs) const
  {
    assert(Rows() == rhs.Rows());
    DenseVector<TYPE> res(Rows());
    for (unsigned int i = 0; i < Rows(); ++i)
      res(i) = (*this)(i) + rhs(i);
    return res;
  }

  /// Subtract matrix from this matrix
  DenseVector<TYPE> operator-(const DenseVector<TYPE>& rhs) const
  {
    assert(Rows() == rhs.Rows());
    DenseVector<TYPE> res(Rows());
    for (unsigned int i = 0; i < Rows(); ++i)
      res(i) = (*this)(i)-rhs(i);
    return res;
  }

  /// Prints the vector to a string and then returns the string.
  std::string PrintStr() const
  {
    std::stringstream out;
    out << "[";
    for (int i = 0; i < (Rows() - 1); ++i)
      out << (*this)(i) << " ";
    out << (*this)(Rows() - 1) << "]";

    return out.str();
  }

  std::vector<TYPE> ToStdVector() const
  {
    std::vector<TYPE> res(Rows());
    for (unsigned int i = 0; i < Rows(); ++i)
      res[i] = (*this)(i);
    return res;
  }
};

/// Multiplication by a scalar from the left.
template <class TYPE>
DenseVector<TYPE>
operator*(const double value, const DenseVector<TYPE>& that)
{
  DenseVector<TYPE> res(that.Rows());
  for (int i = 0; i < that.Rows(); ++i)
    res(i) = that(i) * value;
  return res;
}

} // namespace opensn
