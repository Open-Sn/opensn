#pragma once

#include "framework/mesh/mesh_vector.h"

#include <iostream>

#include <vector>
#include <array>
#include <cmath>
#include <sstream>

#include <type_traits>

namespace opensn
{
struct TensorRank2Dim3;
struct Vector3;
template <int R, int N, class NumberFormat>
struct TensorRNX;

/**Generalized vector notion.
 * \author Jerry, Jan.*/
template <int N, class NumberFormat>
struct VectorNX
{
  std::array<NumberFormat, N> elements;
  const unsigned int dimension;

  /**Default constructor. */
  VectorNX() : dimension(N)
  {
    static_assert(std::is_floating_point<NumberFormat>::value,
                  "Only floating point number formats are supported for VectorNX.");

    elements.fill(NumberFormat());
  }

  /** Constructor with value. */
  VectorNX(const NumberFormat value) : dimension(N)
  {
    static_assert(std::is_floating_point<NumberFormat>::value,
                  "Only floating point number formats are supported for VectorNX.");

    elements.fill(value);
  }

  /** Constructor with Vector3. */
  VectorNX(const Vector3& that) : dimension(N)
  {
    static_assert(std::is_floating_point<NumberFormat>::value,
                  "Only floating point number formats are supported for VectorNX.");

    elements.fill(NumberFormat());

    for (int i = 0; (i < N) and (i < 3); ++i)
      elements[i] = that[i];
  }

  /** Constructor with array of values. This allows constructors of the
   * form:
   * \code
   * Vector3 ihat({1.0,0.0});
   * Vector3 jhat({0.0,1.0,0.0,1.0});
   * Vector3 khat({0.0,0.0,1.0});
   * \endcode */
  VectorNX(const std::vector<NumberFormat>& values) : dimension(N)
  {
    static_assert(std::is_floating_point<NumberFormat>::value,
                  "Only floating point number formats are supported for VectorNX.");

    elements.fill(NumberFormat());
    for (int i = 0; i < values.size() and i < dimension; ++i)
      elements[i] = values[i];
  }

  /**Component-wise copy.*/
  VectorNX& operator=(const VectorNX& rhs)
  {
    elements = rhs.elements;

    return *this;
  }

  /**Component-wise from Vector3.*/
  VectorNX& operator=(const Vector3& rhs)
  {
    for (int i = 0; (i < N) and (i < 3); ++i)
      elements[i] = rhs[i];

    return *this;
  }

  /**Component-wise addition of two vectors.
   * \f$ \vec{w} = \vec{x} + \vec{y} \f$*/
  VectorNX operator+(const VectorNX& rhs) const
  {
    VectorNX<N, NumberFormat> newVector;
    for (int i = 0; i < N; ++i)
      newVector.elements[i] = elements[i] + rhs.elements[i];

    return newVector;
  }

  /**In-place component-wise addition of two vectors.
   * \f$ \vec{x} = \vec{x} + \vec{y} \f$*/
  VectorNX& operator+=(const VectorNX& rhs)
  {
    for (int i = 0; i < N; ++i)
      elements[i] += rhs.elements[i];

    return *this;
  }

  /**Component-wise shift by scalar-value.
   * \f$ \vec{w} = \vec{x} + \alpha \f$*/
  VectorNX Shifted(const NumberFormat value) const
  {
    VectorNX<N, NumberFormat> newVector;
    for (int i = 0; i < N; ++i)
      newVector.elements[i] = elements[i] + value;

    return newVector;
  }

  /**In-place component-wise shift by scalar-value.
   * \f$ \vec{x} = \vec{x} + \alpha \f$*/
  VectorNX& Shift(const NumberFormat value)
  {
    for (int i = 0; i < N; ++i)
      elements[i] += value;

    return *this;
  }

  /**Component-wise subtraction.
   * \f$ \vec{w} = \vec{x} - \vec{y} \f$*/
  VectorNX operator-(const VectorNX& rhs) const
  {
    VectorNX<N, NumberFormat> newVector;
    for (int i = 0; i < N; ++i)
      newVector.elements[i] = elements[i] - rhs.elements[i];

    return newVector;
  }

  /**In-place component-wise subtraction.
   * \f$ \vec{x} = \vec{x} - \vec{y} \f$*/
  VectorNX& operator-=(const VectorNX& rhs)
  {
    for (int i = 0; i < N; ++i)
      elements[i] -= rhs.elements[i];

    return *this;
  }

  /**Vector component-wise multiplication by scalar.
   * \f$ \vec{w} = \vec{x} \alpha \f$*/
  VectorNX operator*(const NumberFormat value) const
  {
    VectorNX<N, NumberFormat> newVector;
    for (int i = 0; i < N; ++i)
      newVector.elements[i] = elements[i] * value;

    return newVector;
  }

  /**Vector in-place component-wise multiplication by scalar.
   * \f$ \vec{x} = \vec{x} \alpha \f$*/
  VectorNX& operator*=(const NumberFormat value)
  {
    for (int i = 0; i < N; ++i)
      elements[i] *= value;

    return *this;
  }

  /**Vector component-wise multiplication.
   * \f$ w_i = x_i y_i \f$*/
  VectorNX operator*(const VectorNX& rhs) const
  {
    VectorNX<N, NumberFormat> newVector;
    for (int i = 0; i < N; ++i)
      newVector.elements[i] = elements[i] * rhs.elements[i];

    return newVector;
  }

  /**Vector in-place component-wise multiplication.
   * \f$ x_i = x_i y_i \f$*/
  VectorNX& operator*=(const VectorNX& rhs)
  {
    for (int i = 0; i < N; ++i)
      elements[i] *= rhs.elements[i];

    return *this;
  }

  /**Vector component-wise division by scalar.
   * \f$ w_i = \frac{x_i}{\alpha} \f$*/
  VectorNX operator/(const NumberFormat value) const
  {
    VectorNX<N, NumberFormat> newVector;
    for (int i = 0; i < N; ++i)
      newVector.elements[i] = elements[i] / value;

    return newVector;
  }

  /**Vector in-place component-wise division by scalar.
   * \f$ x_i = \frac{x_i}{\alpha} \f$*/
  VectorNX& operator/=(const NumberFormat value)
  {
    for (int i = 0; i < N; ++i)
      elements[i] /= value;

    return *this;
  }

  /**Vector component-wise division.
   * \f$ w_i = \frac{x_i}{y_i} \f$*/
  VectorNX operator/(const VectorNX& rhs) const
  {
    VectorNX<N, NumberFormat> newVector;
    for (int i = 0; i < N; ++i)
      newVector.elements[i] = elements[i] / rhs.elements[i];

    return newVector;
  }

  /**Vector in-place component-wise division.
   * \f$ x_i = \frac{x_i}{y_i} \f$*/
  VectorNX& operator/=(const VectorNX& rhs)
  {
    for (int i = 0; i < N; ++i)
      elements[i] /= rhs.elements[i];

    return *this;
  }

  /**Returns a copy of the value at the given index.*/
  NumberFormat operator[](int i) const { return elements[i]; }

  /**Returns a reference of the value at the given index.*/
  NumberFormat& operator()(int i) { return elements[i]; }

  // Tensor product
  TensorRNX<2, N, NumberFormat> OTimes(const VectorNX<N, NumberFormat>& that) const
  {
    TensorRNX<2, N, NumberFormat> out_tensor;
    for (int i = 0; i < N; ++i)
      for (int j = 0; j < N; ++j)
        out_tensor(i)(j) = elements[i] * that.elements[j];

    return out_tensor;
  }

  // Tensor dot product
  VectorNX<N, NumberFormat> Dot(const TensorRNX<2, N, NumberFormat>& that) const
  {
    VectorNX<N, NumberFormat> new_vec;
    for (int i = 0; i < N; ++i)
      new_vec(i) = this->Dot(that.entries[i]);

    return new_vec;
  }

  /**Vector cross-product.
   * \f$ \vec{w} = \vec{x} \times \vec{y} \f$*/
  VectorNX<3, NumberFormat> Cross(const VectorNX<2, NumberFormat>& rhs)
  {
    static_assert(N == 2 or N == 3, "VectorNX::Cross only defined for dimension 2 or 3 vectors.");

    VectorNX<3, NumberFormat> newVector;
    if (dimension == 3)
    {
      newVector(0) = -elements[2] * rhs.elements[1];
      newVector(1) = elements[2] * rhs.elements[0];
      newVector(2) = elements[0] * rhs.elements[1] - elements[1] * rhs.elements[0];
    }
    else
      newVector(2) = elements[0] * rhs.elements[1] - elements[1] * rhs.elements[0];

    return newVector;
  }

  /**Vector cross-product.
   * \f$ \vec{w} = \vec{x} \times \vec{y} \f$*/
  VectorNX<3, NumberFormat> Cross(const VectorNX<3, NumberFormat>& rhs)
  {
    static_assert(N == 2 or N == 3, "VectorNX::Cross only defined for dimension 2 or 3 vectors.");

    VectorNX<3, NumberFormat> newVector;
    if (dimension == 3)
    {
      newVector(0) = elements[1] * rhs.elements[2] - elements[2] * rhs.elements[1];
      newVector(1) = elements[2] * rhs.elements[0] - elements[0] * rhs.elements[2];
      newVector(2) = elements[0] * rhs.elements[1] - elements[1] * rhs.elements[0];
    }
    else
    {
      newVector(0) = elements[1] * rhs.elements[2];
      newVector(1) = -elements[0] * rhs.elements[2];
      newVector(2) = elements[0] * rhs.elements[1] - elements[1] * rhs.elements[0];
    }

    return newVector;
  }

  /**Vector dot-product.
   * \f$ \vec{w} = \vec{x} \bullet \vec{y} \f$ */
  NumberFormat Dot(const VectorNX& rhs) const
  {
    NumberFormat value = 0.0;
    for (int i = 0; i < N; ++i)
      value += elements[i] * rhs.elements[i];

    return value;
  }

  /**Vector dot-product.
   * \f$ \vec{w} = \vec{x} \bullet \vec{y} \f$ */
  NumberFormat Dot(const Vector3& rhs) const
  {
    NumberFormat value = 0.0;
    for (int i = 0; (i < N) and (i < 3); ++i)
      value += elements[i] * rhs[i];

    return value;
  }

  /**Computes the L2-norm of the vector. Otherwise known as the length of
   * a 3D vector.*/
  NumberFormat Norm() const
  {
    NumberFormat value = 0.0;
    for (int i = 0; i < N; ++i)
      value += elements[i] * elements[i];

    value = sqrt(value);
    return value;
  }

  /**Computes the square of the L2-norm of the vector. This eliminates the
   * usage of the square root and is therefore less expensive that a proper
   * L2-norm. Useful if only comparing distances.*/
  NumberFormat NormSquare() const
  {
    NumberFormat value = 0.0;
    for (int i = 0; i < N; ++i)
      value += elements[i] * elements[i];

    return value;
  }

  /**Normalizes the vector in-place.*/
  void Normalize()
  {
    NumberFormat norm = this->Norm();
    for (int i = 0; i < N; ++i)
      elements[i] = elements[i] / norm;
  }

  /**Returns a normalized version of the vector.*/
  VectorNX Normalized() const
  {
    NumberFormat norm = this->Norm();
    VectorNX<N, NumberFormat> newVector;
    for (int i = 0; i < N; ++i)
      newVector.elements[i] = elements[i] / norm;

    return newVector;
  }

  /**Returns a vector v^* where each element is inverted provided
   * that it is greater than the given tolerance, otherwise the offending entry
   * is set to 0.0.*/
  VectorNX InverseZeroIfSmaller(NumberFormat tol) const
  {
    VectorNX<N, NumberFormat> newVector;
    for (int i = 0; i < N; ++i)
      newVector.elements[i] = (std::fabs(elements[i]) > tol) ? 1.0 / elements[i] : 0.0;

    return newVector;
  }

  /**Returns a vector v^* where each element is inverted provided
   * that it is greater than the given tolerance, otherwise the offending entry
   * is set to 1.0.*/
  VectorNX InverseOneIfSmaller(NumberFormat tol) const
  {
    VectorNX<N, NumberFormat> newVector;
    for (int i = 0; i < N; ++i)
      newVector.elements[i] = (std::fabs(elements[i]) > tol) ? 1.0 / elements[i] : 1.0;

    return newVector;
  }

  /**Returns a vector v^* where each element is inverted provided
   * that the inversion is not infinite, otherwise it is zeroed.*/
  VectorNX InverseZeroIfInf() const
  {
    VectorNX<N, NumberFormat> newVector;
    for (int i = 0; i < N; ++i)
    {
      NumberFormat dn_inv = 1.0 / elements[i];
      newVector.elements[i] = (std::isinf(dn_inv)) ? dn_inv : 0.0;
    }

    return newVector;
  }

  /**Returns a vector v^* where each element is inverted without any
   * check for division by zero.
   * \f$ w_i = \frac{1.0}{x_i} \f$*/
  VectorNX Inverse() const
  {
    VectorNX newVector;
    for (int i = 0; i < N; ++i)
      newVector.elements[i] = 1.0 / elements[i];

    return newVector;
  }

  /**prints the vector to standard cout*/
  void Print() const
  {
    for (int i = 0; i < N - 1; i++)
      std::cout << elements[i] << " ";
    std::cout << elements[N - 1];
  }

  // overloading <<

  /**Prints the vector to a string and then returns the string.*/
  std::string PrintS() const
  {
    std::stringstream out;
    out << "[";
    for (int i = 0; i < N - 1; ++i)
      out << elements[i] << " ";
    out << elements[N - 1] << "]";

    return out.str();
  }
};
template <int N>
using VectorN = VectorNX<N, double>;

using Vector2 = VectorN<2>;
// using Vector3 = VectorN<3>;

/**Multiplication by a scalar from the left.*/
template <int N, class NumberFormat>
VectorNX<N, NumberFormat>
operator*(const double value, const VectorNX<N, NumberFormat>& that)
{
  VectorNX<N, NumberFormat> newVector;
  for (int i = 0; i < N; ++i)
    newVector.elements[i] = that.elements[i] * value;
  return newVector;
}

} // namespace opensn
