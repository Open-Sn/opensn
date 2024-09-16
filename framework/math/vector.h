// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include "framework/data_types/ndarray.h"
#include <sstream>
#include <cmath>

namespace opensn
{

template <typename TYPE>
class Vector : public NDArray<TYPE>
{
public:
  /// Create an empty dense (column) vector
  Vector() : NDArray<TYPE>() {}

  /// Create a dense (column) vector with specified number of rows
  Vector(unsigned int rows) : NDArray<TYPE>({rows}) {}

  /// Create a dense (column) vector with specified number of rows and initialize the elements
  Vector(unsigned int rows, TYPE intitial_value) : NDArray<TYPE>({rows}, intitial_value) {}

  Vector(const std::vector<TYPE>& in) : NDArray<TYPE>({in.size()})
  {
    for (auto i = 0; i < in.size(); ++i)
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
      for (auto i = old_size; i < Rows(); ++i)
        (*this)(i) = value;
    }
  }

  /// Set the elements of the vector to a specified value
  void Set(TYPE val) { this->set(val); }

  /// Scale the vector with a constant value
  void Scale(TYPE alpha)
  {
    for (auto i = 0; i < Rows(); ++i)
      (*this)(i) *= alpha;
  }

  /// Normalizes the vector in-place.
  void Normalize()
  {
    TYPE mag = this->Magnitude();
    for (auto i = 0; i < Rows(); ++i)
      (*this)(i) /= mag;
  }

  /// Computes the L2-norm of the vector. Otherwise known as the length of a 3D vector.
  TYPE Magnitude() const
  {
    TYPE value = 0.0;
    for (auto i = 0; i < Rows(); ++i)
      value += (*this)(i) * (*this)(i);
    value = sqrt(value);
    return value;
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
    for (auto i = 0; i < Rows(); ++i)
      res[i] = (*this)(i);
    return res;
  }
};

/// Scale the vector with a constant value
template <typename TYPE>
void
Scale(Vector<TYPE>& a, TYPE alpha)
{
  for (auto i = 0; i < a.Rows(); ++i)
    a(i) *= alpha;
}

template <typename TYPE>
Vector<TYPE>
Scaled(const Vector<TYPE>& a, TYPE alpha)
{
  Vector<TYPE> res(a.Rows());
  for (auto i = 0; i < a.Rows(); ++i)
    res(i) = a(i) * alpha;
  return res;
}

/// Add vector to this vector
template <typename TYPE>
Vector<TYPE>
Add(const Vector<TYPE>& a, const Vector<TYPE>& b)
{
  assert(a.Rows() == b.Rows());
  Vector<TYPE> res(a.Rows());
  for (auto i = 0; i < a.Rows(); ++i)
    res(i) = a(i) + b(i);
  return res;
}

/// Subtract two vectors
template <typename TYPE>
Vector<TYPE>
Subtract(const Vector<TYPE>& a, const Vector<TYPE>& b)
{
  assert(a.Rows() == b.Rows());
  Vector<TYPE> res(a.Rows());
  for (auto i = 0; i < a.Rows(); ++i)
    res(i) = a(i) - b(i);
  return res;
}

template <typename TYPE>
double
Vec2Norm(const Vector<TYPE>& x)
{
  auto n = x.Rows();
  double val = 0.0;
  for (auto i = 0; i != n; ++i)
    val += x(i) * x(i);
  return std::sqrt(val);
}

template <typename TYPE>
double
Dot(const Vector<TYPE>& x, const Vector<TYPE>& y)
{
  assert(x.Rows() > 0);
  assert(y.Rows() > 0);
  assert(x.Rows() == y.Rows());
  double val = 0.0;
  for (auto i = 0; i < x.Rows(); ++i)
    val += x(i) * y(i);
  return val;
}

} // namespace opensn
