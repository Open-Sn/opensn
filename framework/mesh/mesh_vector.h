// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include <vector>
#include <array>
#include <iostream>
#include <cmath>
#include <sstream>
#include <memory>

namespace opensn
{

struct TensorRank2Dim3;

/// General 3-element vector structure.
struct Vector3 : public std::enable_shared_from_this<Vector3>
{
  /// X-component of the vector
  double x{0.0};
  /// Y-component of the vector
  double y{0.0};
  /// Z-component of the vector
  double z{0.0};

  /// Default constructor: initializes to (0, 0, 0).
  Vector3() = default;

  /// Constructor where \f$ \vec{x}=[a,b,c] \f$.
  explicit Vector3(double a, double b = 0.0, double c = 0.0) : x(a), y(b), z(c) {}

  /// Constructor where \f$ \vec{x}=\{a,b,c\} \f$.
  Vector3(std::initializer_list<double> list)
  {
    auto it = list.begin();
    if (list.size() > 0)
      x = *it++;
    if (list.size() > 1)
      y = *it++;
    if (list.size() > 2)
      z = *it;
  }

  /// Constructor where \f$ \vec{x}=\{a,b,c\} \f$.
  explicit Vector3(const std::vector<double>& list)
  {
    if (!list.empty())
      x = list[0];
    if (list.size() > 1)
      y = list[1];
    if (list.size() > 2)
      z = list[2];
  }

  /// Component-wise addition of two vectors. \f$ \vec{w} = \vec{x} + \vec{y} \f$
  Vector3 operator+(const Vector3& other) const { return {x + other.x, y + other.y, z + other.z}; }

  /// In-place component-wise addition of two vectors. \f$ \vec{x} = \vec{x} + \vec{y} \f$
  Vector3& operator+=(const Vector3& that)
  {
    x += that.x;
    y += that.y;
    z += that.z;
    return *this;
  }

  /// Component-wise shift by scalar-value. \f$ \vec{w} = \vec{x} + \alpha \f$
  Vector3 Shifted(double value) const { return Vector3{x + value, y + value, z + value}; }

  /// In-place component-wise shift by scalar value. \f$ \vec{x} = \vec{x} + \alpha \f$
  Vector3& Shift(double value)
  {
    x += value;
    y += value;
    z += value;
    return *this;
  }

  /// Component-wise subtraction of two vectors. \f$ \vec{w} = \vec{x} - \vec{y} \f$
  Vector3 operator-(const Vector3& other) const { return {x - other.x, y - other.y, z - other.z}; }

  /// In-place component-wise subtraction. \f$ \vec{x} = \vec{x} - \vec{y} \f$
  Vector3& operator-=(const Vector3& that)
  {
    x -= that.x;
    y -= that.y;
    z -= that.z;
    return *this;
  }

  /// Vector component-wise multiplication by scalar. \f$ \vec{w} = \vec{x} \alpha \f$
  Vector3 operator*(double value) const { return Vector3{x * value, y * value, z * value}; }

  /// Vector in-place component-wise multiplication by scalar. \f$ \vec{x} = \vec{x} \alpha \f$
  Vector3& operator*=(double value)
  {
    x *= value;
    y *= value;
    z *= value;
    return *this;
  }

  /// Vector component-wise multiplication. \f$ w_i = x_i y_i \f$
  Vector3 operator*(const Vector3& that) const
  {
    return Vector3{x * that.x, y * that.y, z * that.z};
  }

  /// Vector in-place component-wise multiplication. \f$ x_i = x_i y_i \f$
  Vector3& operator*=(const Vector3& that)
  {
    x *= that.x;
    y *= that.y;
    z *= that.z;
    return *this;
  }

  /// Vector component-wise division by scalar. \f$ w_i = \frac{x_i}{\alpha} \f$
  Vector3 operator/(double value) const
  {
    if (value == 0.0)
      throw std::runtime_error("Division by zero in Vector3 operator/");
    return Vector3{x / value, y / value, z / value};
  }

  /// Vector in-place component-wise division by scalar. \f$ x_i = \frac{x_i}{\alpha} \f$
  Vector3& operator/=(double value)
  {
    if (value == 0.0)
      throw std::runtime_error("Division by zero in Vector3 operator/=");
    x /= value;
    y /= value;
    z /= value;
    return *this;
  }

  /// Vector component-wise division. \f$ w_i = \frac{x_i}{y_i} \f$
  Vector3 operator/(const Vector3& that) const
  {
    if (that.x == 0.0 || that.y == 0.0 || that.z == 0.0)
      throw std::runtime_error("Division by zero in Vector3 operator/ (component-wise)");
    return Vector3{x / that.x, y / that.y, z / that.z};
  }

  /// Vector in-place component-wise division. \f$ x_i = \frac{x_i}{y_i} \f$
  Vector3& operator/=(const Vector3& that)
  {
    if (that.x == 0.0 || that.y == 0.0 || that.z == 0.0)
      throw std::runtime_error("Division by zero in Vector3 operator/= (component-wise)");
    x /= that.x;
    y /= that.y;
    z /= that.z;
    return *this;
  }

  /// Returns a copy of the value at the given index.
  double operator[](size_t i) const
  {
    if (i > 2)
      throw std::out_of_range("Index out of range in Vector3 operator[]");
    return i == 0 ? x : (i == 1 ? y : z);
  }

  /// Returns a reference of the value at the given index.
  double& operator()(size_t i)
  {
    if (i > 2)
      throw std::out_of_range("Index out of range in Vector3 operator()");
    return i == 0 ? x : (i == 1 ? y : z);
  }

  /**
   * Vector cross-product.
   * \f$ \vec{w} = \vec{x} \times \vec{y} \f$
   */
  Vector3 Cross(const Vector3& that) const
  {
    Vector3 newVector;
    newVector.x = this->y * that.z - this->z * that.y;
    newVector.y = this->z * that.x - this->x * that.z;
    newVector.z = this->x * that.y - this->y * that.x;

    return newVector;
  }

  /// Vector dot-product. \f$ \vec{w} = \vec{x} \bullet \vec{y} \f$
  double Dot(const Vector3& that) const { return x * that.x + y * that.y + z * that.z; }

  /// Computes the L2-norm of the vector. Otherwise known as the length of a 3D vector.
  double Norm() const { return std::sqrt(NormSquare()); }

  /// Computes the square of the L2-norm of the vector. Less expensive than a proper L2-norm.
  double NormSquare() const { return x * x + y * y + z * z; }

  /// Normalizes the vector in-place. \f$ \vec{x} = \frac{\vec{x}}{||x||_2} \f$
  void Normalize()
  {
    double norm = Norm();
    if (norm == 0.0)
      throw std::runtime_error("Cannot normalize a zero-length vector.");
    x /= norm;
    y /= norm;
    z /= norm;
  }

  /// Returns a normalized version of the vector. \f$ \vec{w} = \frac{\vec{x}}{||x||_2} \f$
  Vector3 Normalized() const
  {
    double norm = Norm();
    if (norm == 0.0)
      throw std::runtime_error("Cannot normalize a zero-length vector.");
    return Vector3{x / norm, y / norm, z / norm};
  }

  /**
   * Returns a vector v^* where each element is inverted provided that it is greater than the given
   * tolerance, otherwise the offending entry is set to 0.0. \f$ w_i = \frac{1.0}{x_i} \f$
   */
  /// Inverse components, setting zero if below tolerance.
  Vector3 InverseZeroIfSmaller(double tol) const
  {
    return Vector3{std::fabs(x) > tol ? 1.0 / x : 0.0,
                   std::fabs(y) > tol ? 1.0 / y : 0.0,
                   std::fabs(z) > tol ? 1.0 / z : 0.0};
  }

  /// Inverse components, setting one if below tolerance.
  Vector3 InverseOneIfSmaller(double tol) const
  {
    return Vector3{std::fabs(x) > tol ? 1.0 / x : 1.0,
                   std::fabs(y) > tol ? 1.0 / y : 1.0,
                   std::fabs(z) > tol ? 1.0 / z : 1.0};
  }

  /// Inverts each component without checking for division by zero.
  Vector3 Inverse() const
  {
    if (x == 0.0 || y == 0.0 || z == 0.0)
      throw std::runtime_error("Division by zero in Vector3::Inverse.");
    return Vector3{1.0 / x, 1.0 / y, 1.0 / z};
  }

  /// Prints the vector to std::cout.
  void Print() const { std::cout << "[" << x << " " << y << " " << z << "]"; }

  /// Returns the vector as a string.
  std::string PrintStr() const
  {
    std::ostringstream out;
    out << "[" << x << " " << y << " " << z << "]";
    return out.str();
  }

  static constexpr size_t Size() { return 3; }
};

/// Scalar multiplication from the left. \f$ \vec{w} = \alpha \vec{x}\f$
Vector3 operator*(double value, const Vector3& that);

} // namespace opensn
