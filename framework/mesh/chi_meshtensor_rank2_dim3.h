#pragma once

#include "framework/mesh/chi_meshvector.h"

/**General rank 2 tensor to be used with Vector3*/
struct chi_mesh::TensorRank2Dim3
{
  std::vector<chi_mesh::Vector3> t; ///< Tensor entries.

  /**Default constructor.*/
  TensorRank2Dim3() { t.resize(3, chi_mesh::Vector3(0.0, 0.0, 0.0)); }
  /**Constructor with specified value.
   * \f$ T_{ii} = \alpha \f$*/
  explicit TensorRank2Dim3(const double value)
  {
    t.resize(3, chi_mesh::Vector3(value, value, value));
  }

  /**Copy constructor.*/
  TensorRank2Dim3(const TensorRank2Dim3& that) { this->t = that.t; }

  /**Component-wise copy*/
  TensorRank2Dim3& operator=(const TensorRank2Dim3& that)
  {
    this->t = that.t;

    return *this;
  }

  /**Element access.*/
  chi_mesh::Vector3& operator[](int index) { return t[index]; }

  /**Component-wise addition.
   * \f$ \vec{\vec{W}} = \vec{\vec{X}} + \vec{\vec{Y}} \f$*/
  TensorRank2Dim3 operator+(const TensorRank2Dim3& that) const
  {
    TensorRank2Dim3 new_t;
    for (int i = 0; i < 3; ++i)
      for (int j = 0; j < 3; ++j)
        new_t[i](j) = this->t[i][j] + that.t[i][j];

    return new_t;
  }

  /**In-place component-wise addition.
   * \f$ \vec{\vec{X}} = \vec{\vec{X}} + \vec{\vec{Y}} \f$*/
  TensorRank2Dim3& operator+=(const TensorRank2Dim3& that)
  {
    for (int i = 0; i < 3; ++i)
      for (int j = 0; j < 3; ++j)
        this->t[i](j) += that.t[i][j];

    return *this;
  }

  /**Component-wise subtraction.
   * \f$ \vec{\vec{W}} = \vec{\vec{X}} - \vec{\vec{Y}} \f$*/
  TensorRank2Dim3 operator-(const TensorRank2Dim3& that) const
  {
    TensorRank2Dim3 new_t;
    for (int i = 0; i < 3; ++i)
      for (int j = 0; j < 3; ++j)
        new_t[i](j) = this->t[i][j] - that.t[i][j];

    return new_t;
  }

  /**In-place component-wise subtraction.
   * \f$ \vec{\vec{X}} = \vec{\vec{X}} - \vec{\vec{Y}} \f$*/
  TensorRank2Dim3& operator-=(const TensorRank2Dim3& that)
  {
    for (int i = 0; i < 3; ++i)
      for (int j = 0; j < 3; ++j)
        this->t[i](j) -= that.t[i][j];

    return *this;
  }

  /**Component-wise multiplication by scalar.
   * \f$ \vec{\vec{W}} = \vec{\vec{X}}\alpha \f$ */
  TensorRank2Dim3 operator*(const double value) const
  {
    TensorRank2Dim3 new_t;
    for (int i = 0; i < 3; ++i)
      for (int j = 0; j < 3; ++j)
        new_t[i](j) = this->t[i][j] * value;

    return new_t;
  }

  /**In-place component-wise multiplication by scalar.
   *  \f$ \vec{\vec{X}} = \vec{\vec{X}}\alpha \f$*/
  TensorRank2Dim3& operator*=(const double value)
  {
    for (int i = 0; i < 3; ++i)
      for (int j = 0; j < 3; ++j)
        this->t[i](j) *= value;

    return *this;
  }

  /**Component-wise division by scalar.
   *  \f$ \vec{\vec{W}} = \vec{\vec{X}}\frac{1}{\alpha} \f$*/
  TensorRank2Dim3 operator/(const double value) const
  {
    TensorRank2Dim3 new_t;
    for (int i = 0; i < 3; ++i)
      for (int j = 0; j < 3; ++j)
        new_t[i](j) = this->t[i][j] / value;

    return new_t;
  }

  /**In-place component-wise division by scalar.
   * \f$ \vec{\vec{X}} = \vec{\vec{X}}\frac{1}{\alpha} \f$*/
  TensorRank2Dim3& operator/=(const double value)
  {
    for (int i = 0; i < 3; ++i)
      for (int j = 0; j < 3; ++j)
        this->t[i](j) /= value;

    return *this;
  }

  /**Classical transpose of the tensor.
   * \f$ W_{ij} = T_{ji} \f$*/
  TensorRank2Dim3 Transpose()
  {
    TensorRank2Dim3 new_t;
    for (int i = 0; i < 3; ++i)
      for (int j = 0; j < 3; ++j)
        new_t[i](j) = this->t[j][i];

    return new_t;
  }

  /**
   * Dot product of rank-2 tensor with a vector.
   * \f$ \vec{w} = \vec{\vec{T}} \bullet \vec{x} \f$
   */
  Vector3 Dot(const chi_mesh::Vector3& v) const;

  /**
   * Returns the diagonal of a rank-2 dim-3 tensor as a vector3.
   * \f$ \vec{w} = \text{diag} \vec{\vec{T}} \f$
   */
  chi_mesh::Vector3 Diag() const;

  /**Returns the sum of the diagonal. Sometimes useful to get
   * divergence of a vector given its gradient.*/
  double DiagSum() const
  {
    double val = 0.0;
    for (int i = 0; i < 3; ++i)
      val += t[i][i];

    return val;
  }

  /**Prints the vector to a string and then returns the string.*/
  std::string PrintS()
  {
    std::stringstream out;
    out << "[";
    for (int i = 0; i < 3; i++)
    {
      for (int j = 0; j < 3; j++)
      {
        out << t[i][j] << " ";
      }
      if (i != 2) out << "\n ";
      else
        out << "]";
    }

    return out.str();
  }
};

/**
 * Rank-2 dim-3 tensor multiplied from the left with a scalar.
 * \f$ \vec{\vec{W}} = \alpha \vec{\vec{T}} \f$
 */
chi_mesh::TensorRank2Dim3 operator*(const double value, const chi_mesh::TensorRank2Dim3& that);
