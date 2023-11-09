#pragma once

#include "framework/math/math_vector_nx.h"

#include <iostream>

#include <vector>
#include <cmath>
#include <sstream>

#include <type_traits>

namespace chi_math
{

/**Generalized tensor with rank-R dimension-N and arbitrary number
 * format.
 * \author Jan.*/
template <int R, int N, class NumberFormat>
struct TensorRNX
{
  std::vector<VectorNX<N, NumberFormat>> entries;
  const unsigned int rank;
  const unsigned int dimension;

  /**Default constructor.*/
  TensorRNX() : rank(R), dimension(N)
  {
    static_assert(std::is_floating_point<NumberFormat>::value,
                  "Only floating point number formats are supported for TensorRNX.");

    const int num_rank1_entries = std::pow(N, R - 1);
    entries.resize(num_rank1_entries);
  }

  /**Constructor with value.*/
  TensorRNX(const NumberFormat value) : rank(R), dimension(N)
  {
    static_assert(std::is_floating_point<NumberFormat>::value,
                  "Only floating point number formats are supported for TensorRNX.");

    const int num_rank1_entries = std::pow(N, R - 1);
    entries.resize(num_rank1_entries, VectorNX<N, NumberFormat>(value));
  }

  /**Component-wise copy.*/
  TensorRNX& operator=(const TensorRNX& rhs)
  {
    entries = rhs.entries;

    return *this;
  }

  template <int R2>
  struct RecursiveAccessor
  {
    const int int_rank;

    RecursiveAccessor(int i) : int_rank(R2) {}
  };

  VectorNX<N, NumberFormat> operator[](const int i);

  /***/
  //    RecursiveAccessor operator[](const int i)
  //    {
  //
  //    }
};

/**Specialized rank-2 tensor.*/
template <int N, class NumberFormat>
struct TensorRNX<2, N, NumberFormat>
{
  std::vector<VectorNX<N, NumberFormat>> entries;
  const unsigned int rank;
  const unsigned int dimension;

  /**Default constructor.*/
  TensorRNX() : rank(2), dimension(N)
  {
    static_assert(std::is_floating_point<NumberFormat>::value,
                  "Only floating point number formats are supported for TensorRNX.");

    const int num_rank1_entries = std::pow(N, 2 - 1);
    entries.resize(num_rank1_entries);
  }

  /**Constructor with value.*/
  TensorRNX(const NumberFormat value) : rank(2), dimension(N)
  {
    static_assert(std::is_floating_point<NumberFormat>::value,
                  "Only floating point number formats are supported for TensorRNX.");

    const int num_rank1_entries = std::pow(N, 2 - 1);
    entries.resize(num_rank1_entries, VectorNX<N, NumberFormat>(value));
  }

  /**Copy constructor.*/
  TensorRNX(const TensorRNX<2, N, NumberFormat>& that) : rank(2), dimension(N)
  {
    this->entries = that.entries;
  }

  /**Component-wise copy (assignment operator.*/
  TensorRNX& operator=(const TensorRNX& rhs)
  {
    entries = rhs.entries;

    return *this;
  }

  /**Return reference to vector at given row.*/
  const VectorNX<N, NumberFormat>& operator[](const int i) const { return entries.at(i); }

  /**Return reference to vector at given row.*/
  VectorNX<N, NumberFormat>& operator()(const int i) { return entries.at(i); }

  /**Component-wise addition.
   * \f$ \vec{\vec{W}} = \vec{\vec{X}} + \vec{\vec{Y}} \f$*/
  TensorRNX<2, N, NumberFormat> operator+(const TensorRNX<2, N, NumberFormat>& that) const
  {
    TensorRNX<2, N, NumberFormat> new_t;
    for (int i = 0; i < N; ++i)
    {
      const auto& this_row = entries[i];
      const auto& that_row = that.entries[i];
      for (int j = 0; j < N; ++j)
        new_t(i)(j) = this_row[j] + that_row[j];
    }

    return new_t;
  }

  /**In-place component-wise addition.
   * \f$ \vec{\vec{X}} = \vec{\vec{X}} + \vec{\vec{Y}} \f$*/
  TensorRNX<2, N, NumberFormat>& operator+=(const TensorRNX<2, N, NumberFormat>& that)
  {
    for (int i = 0; i < N; ++i)
    {
      const auto& that_row = that.entries[i];
      for (int j = 0; j < N; ++j)
        this->entries[i](j) += that_row[j];
    }

    return *this;
  }

  /**Component-wise subtraction.
   * \f$ \vec{\vec{W}} = \vec{\vec{X}} - \vec{\vec{Y}} \f$*/
  TensorRNX<2, N, NumberFormat> operator-(const TensorRNX<2, N, NumberFormat>& that) const
  {
    TensorRNX<2, N, NumberFormat> new_t;
    for (int i = 0; i < N; ++i)
    {
      const auto& this_row = entries[i];
      const auto& that_row = that.entries[i];
      for (int j = 0; j < N; ++j)
        new_t(i)(j) = this_row[j] - that_row[j];
    }

    return new_t;
  }

  /**In-place component-wise subtraction.
   * \f$ \vec{\vec{X}} = \vec{\vec{X}} - \vec{\vec{Y}} \f$*/
  TensorRNX<2, N, NumberFormat>& operator-=(const TensorRNX<2, N, NumberFormat>& that)
  {
    for (int i = 0; i < N; ++i)
    {
      const auto& that_row = that.entries[i];
      for (int j = 0; j < N; ++j)
        this->entries[i](j) -= that_row[j];
    }

    return *this;
  }

  /**Component-wise multiplication by scalar.
   * \f$ \vec{\vec{W}} = \vec{\vec{X}}\alpha \f$ */
  TensorRNX<2, N, NumberFormat> operator*(const double value) const
  {
    TensorRNX<2, N, NumberFormat> new_t;
    for (int i = 0; i < N; ++i)
    {
      const auto& this_row = entries[i];
      for (int j = 0; j < N; ++j)
        new_t(i)(j) = this_row[j] * value;
    }

    return new_t;
  }

  /**In-place component-wise multiplication by scalar.
   *  \f$ \vec{\vec{X}} = \vec{\vec{X}}\alpha \f$*/
  TensorRNX<2, N, NumberFormat>& operator*=(const double value)
  {
    for (int i = 0; i < N; ++i)
      for (int j = 0; j < N; ++j)
        this->entries[i](j) *= value;

    return *this;
  }

  /**Component-wise division by scalar.
   *  \f$ \vec{\vec{W}} = \vec{\vec{X}}\frac{1}{\alpha} \f$*/
  TensorRNX<2, N, NumberFormat> operator/(const double value) const
  {
    TensorRNX<2, N, NumberFormat> new_t;
    for (int i = 0; i < N; ++i)
    {
      const auto& this_row = entries[i];
      for (int j = 0; j < N; ++j)
        new_t(i)(j) = this_row[j] / value;
    }

    return new_t;
  }

  /**In-place component-wise division by scalar.
   * \f$ \vec{\vec{X}} = \vec{\vec{X}}\frac{1}{\alpha} \f$*/
  TensorRNX<2, N, NumberFormat>& operator/=(const double value)
  {
    for (int i = 0; i < N; ++i)
      for (int j = 0; j < N; ++j)
        this->entries[i](j) /= value;

    return *this;
  }

  /**Classical transpose of the tensor.
   * \f$ W_{ij} = T_{ji} \f$*/
  TensorRNX<2, N, NumberFormat> Transpose() const
  {
    TensorRNX<2, N, NumberFormat> new_t;
    for (int i = 0; i < N; ++i)
      for (int j = 0; j < N; ++j)
      {
        const auto& this_row = entries[j];
        new_t(i)(j) = this_row[i];
      }

    return new_t;
  }

  /**Tensor dot-product with a vector.*/
  VectorNX<N, NumberFormat> Dot(const VectorNX<N, NumberFormat>& v) const
  {
    VectorNX<N, NumberFormat> new_vec;
    for (int i = 0; i < N; ++i)
    {
      const auto& this_row = entries[i];
      new_vec(i) = this_row.Dot(v);
    }

    return new_vec;
  }

  /**Obtains the diagonal of a rank-2 tensor as a vector.*/
  VectorNX<N, NumberFormat> Diag() const
  {
    VectorNX<N, NumberFormat> new_vec;
    for (int i = 0; i < N; ++i)
    {
      const auto& this_row = entries[i];
      new_vec(i) = this_row[i];
    }

    return new_vec;
  }

  /**Returns the sum of the diagonal. Sometimes useful to get
   * divergence of a vector given its gradient.*/
  double DiagSum() const
  {
    double val = 0.0;
    for (int i = 0; i < N; ++i)
    {
      const auto& this_row = entries[i];
      val += this_row[i];
    }

    return val;
  }

  /**Prints the tensor to a string and then returns the string.*/
  std::string PrintS() const
  {
    std::stringstream out;
    out << "[";
    for (int i = 0; i < N; ++i)
    {
      const auto& row = this->operator[](i);
      for (int j = 0; j < (N - 1); ++j)
        out << row[j] << " ";
      if (i == (N - 1)) out << row[N - 1] << "]";
      else
        out << row[N - 1] << "\n";
    }

    return out.str();
  }
};
template <int R, int N>
using TensorRN = TensorRNX<R, N, double>;

template <int N>
using Tensor2N = TensorRNX<2, N, double>;

} // namespace chi_math

/**Multiplication by a scalar from the left.*/
template <int N, class NumberFormat>
chi_math::TensorRNX<2, N, NumberFormat>
operator*(const double value, const chi_math::TensorRNX<2, N, NumberFormat>& that)
{
  chi_math::TensorRNX<2, N, NumberFormat> new_tensor;
  for (int i = 0; i < N; ++i)
    for (int j = 0; j < N; ++j)
      new_tensor(i)(j) = that[i][j] * value;

  return new_tensor;
}
