// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include <cmath>
#include <sstream>
#include <array>
#include <algorithm>
#include "framework/math/vector3.h"

namespace opensn
{

struct Matrix3x3
{
  /// Storage for matrix elements
  std::array<double, 9> vals{0., 0., 0., 0., 0., 0., 0., 0., 0.};

  /**
   * Produces a rotation matrix with a reference vector rotated from the
   * cartesian basis vectors \f$\hat{i}\f$, \f$\hat{j}\f$ and \f$\hat{k}\f$.
   *
   * By default a rotation matrix that creates no rotation is
   * the identity matrix. Such a matrix can be defined from basis vectors
   * following the notion that the "up-vector" is \f$\hat{k}\f$,
   * this is also called the normal vector \f$\hat{n}\f$.
   * The tangent vector is \f$\hat{i}\f$, denoted with \f$\hat{t}\f$.
   * And the bi-norm vector is \f$\hat{j}\f$, denoted with \f$\hat{b}\f$.
   *
   * By specifying only the normal vector we can compute a simple pitch based
   * rotation matrix. The supplied vector is therefore the new normal-vector,
   * the tangent vector is computed as \f$ \hat{t} = \hat{n} \times \hat{k} \f$,
   * and the bi-norm vector is computed as
   * \f$ \hat{b} = \hat{n} \times \hat{t} \f$
   */
  static Matrix3x3 MakeRotationMatrixFromVector(const Vector3& vec)
  {
    Matrix3x3 R;

    Vector3 n = vec;
    Vector3 khat(0.0, 0.0, 1.0);

    if (n.Dot(khat) > 0.9999999)
      R.SetDiagonalVec(1.0, 1.0, 1.0);
    else if (n.Dot(khat) < -0.9999999)
      R.SetDiagonalVec(1.0, 1.0, -1.0);
    else
    {
      auto tangent = n.Cross(khat).Normalized();
      auto binorm = n.Cross(tangent).Normalized();

      R.SetColJVec(0, tangent);
      R.SetColJVec(1, binorm);
      R.SetColJVec(2, n);
    }
    return R;
  }

  /// Matrix addition operator
  Matrix3x3 operator+(const Matrix3x3& other) const
  {
    Matrix3x3 result;
    std::transform(
      vals.begin(), vals.end(), other.vals.begin(), result.vals.begin(), std::plus<>());
    return result;
  }

  /// Matrix subtraction operator
  Matrix3x3 operator-(const Matrix3x3& other) const
  {
    Matrix3x3 result;
    std::transform(
      vals.begin(), vals.end(), other.vals.begin(), result.vals.begin(), std::minus<>());
    return result;
  }

  /// Matrix scalar multiplication operator
  Matrix3x3 operator*(const double scalar) const
  {
    Matrix3x3 result;
    std::transform(
      vals.begin(), vals.end(), result.vals.begin(), [scalar](double val) { return val * scalar; });
    return result;
  }

  /// In-place scalar multiplication operator
  Matrix3x3& operator*=(const double scalar)
  {
    for (auto& val : vals)
      val *= scalar;
    return *this;
  }

  /// Matrix multiply with vector.
  Vector3 operator*(const Vector3& vec) const
  {
    std::array<double, 3> i_vec = {vec.x, vec.y, vec.z};
    std::array<double, 3> o_vec = {0.0, 0.0, 0.0};

    for (int i = 0; i < 3; ++i)
      for (int j = 0; j < 3; ++j)
        o_vec[i] += this->GetIJ(i, j) * i_vec[j];

    return Vector3(o_vec[0], o_vec[1], o_vec[2]);
  }

  /// Set value at row i and column j.
  void SetIJ(int i, int j, double value) { vals[i * 3 + j] = value; }

  /// Add value to value at row i and column j.
  void AddIJ(int i, int j, double value) { vals[i * 3 + j] += value; }

  /// Obtain a copy of the value at row i and column j.
  double GetIJ(int i, int j) const { return vals[i * 3 + j]; }

  /// Set row i using a vector.
  void SetRowIVec(int i, const Vector3& vec)
  {
    vals[i * 3 + 0] = vec.x;
    vals[i * 3 + 1] = vec.y;
    vals[i * 3 + 2] = vec.z;
  }

  /// Set column j using a vector.
  void SetColJVec(int j, const Vector3& vec)
  {
    vals[0 * 3 + j] = vec.x;
    vals[1 * 3 + j] = vec.y;
    vals[2 * 3 + j] = vec.z;
  }

  /// Sets the diagonal of the matrix.
  void SetDiagonalVec(double a00, double a11, double a22)
  {
    vals[0 * 3 + 0] = a00;
    vals[1 * 3 + 1] = a11;
    vals[2 * 3 + 2] = a22;
  }

  /// Get the determinant using specified row [default:0].
  double Det(int row = 0) const
  {
    double det = 0.0;
    int sign = -1;

    for (int j = 0; j < 3; ++j)
    {
      sign *= -1;
      det += sign * GetIJ(row, j) * MinorIJ(row, j);
    }

    return det;
  }

  /// Get the minor value associated with row ir and column jr.
  double MinorIJ(int ir, int jr) const
  {
    std::array<double, 4> a{};
    int k = 0;

    for (int i = 0; i < 3; ++i)
    {
      if (i == ir)
        continue;
      for (int j = 0; j < 3; ++j)
      {
        if (j == jr)
          continue;
        a[k++] = vals[i * 3 + j];
      }
    }

    return a[0] * a[3] - a[1] * a[2];
  }

  /// Compute the matrix transpose.
  Matrix3x3 Transpose() const
  {
    Matrix3x3 result;
    for (int i = 0; i < 3; ++i)
      for (int j = 0; j < 3; ++j)
        result.vals[j * 3 + i] = vals[i * 3 + j]; // Swap row and column

    return result;
  }

  /// Compute the matrix inverse.
  Matrix3x3 Inverse() const
  {
    Matrix3x3 oM;
    Matrix3x3 oMT;

    // Compute matrix of minors
    for (int i = 0; i < 3; ++i)
      for (int j = 0; j < 3; ++j)
        oM.SetIJ(i, j, MinorIJ(i, j));

    // Compute matrix of cofactors
    int sign = -1;
    for (int i = 0; i < 3; ++i)
      for (int j = 0; j < 3; ++j)
      {
        sign *= -1;
        oM.SetIJ(i, j, oM.GetIJ(i, j) * sign);
      }

    // Compute the transpose
    oMT = oM.Transpose();

    // Get determinant
    double det = Det();

    if (det == 0.0)
      throw std::runtime_error("Matrix is singular and cannot be inverted.");

    return oMT * (1.0 / det); // Multiply by reciprocal of determinant
  }

  /// Print the matrix as a string.
  std::string PrintStr() const
  {
    std::stringstream out;
    out << "[";
    for (int i = 0; i < 3; ++i)
    {
      for (int j = 0; j < 3; ++j)
      {
        out << GetIJ(i, j) << " ";
      }
      if (i != 2)
        out << "\n ";
      else
        out << "]";
    }

    return out.str();
  }
};

} // namespace opensn
