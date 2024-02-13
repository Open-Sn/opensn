#pragma once

#include <cmath>
#include <sstream>
#include "framework/mesh/mesh_vector.h"

namespace opensn
{

struct Matrix3x3
{
  double vals[9];

  /**Creates a zeros matrix of size 3x3.*/
  Matrix3x3()
  {
    vals[0] = 0.0;
    vals[1] = 0.0;
    vals[2] = 0.0;
    vals[3] = 0.0;
    vals[4] = 0.0;
    vals[5] = 0.0;
    vals[6] = 0.0;
    vals[7] = 0.0;
    vals[8] = 0.0;
  }

  /**Produces a rotation matrix with a reference vector rotated from the
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
   * \f$ \hat{b} = \hat{n} \times \hat{t} \f$*/
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

  /**Copy constructor*/
  Matrix3x3& operator=(const Matrix3x3& inM)
  {
    for (int k = 0; k < 9; k++)
      this->vals[k] = inM.vals[k];
    return *this;
  }

  /**Matrix addition operator.*/
  Matrix3x3 operator+(const Matrix3x3& inM)
  {
    Matrix3x3 oM;
    for (int k = 0; k < 9; k++)
      oM.vals[k] = this->vals[k] + inM.vals[k];
    return oM;
  }

  /**Matrix subtraction operator.*/
  Matrix3x3 operator-(const Matrix3x3& inM)
  {
    Matrix3x3 oM;
    for (int k = 0; k < 9; k++)
      oM.vals[k] = this->vals[k] - inM.vals[k];
    return oM;
  }

  /**Matrix multiply with scalar.*/
  Matrix3x3 operator*(const double value)
  {
    Matrix3x3 oM;
    for (int k = 0; k < 9; k++)
    {
      oM.vals[k] = this->vals[k] * value;
    }
    return oM;
  }

  /**Matrix multiply with vector.*/
  Vector3 operator*(const Vector3& vec) const
  {
    double i_vec[] = {vec.x, vec.y, vec.z};
    double o_vec[] = {0.0, 0.0, 0.0};

    for (int i = 0; i < 3; i++)
    {
      for (int j = 0; j < 3; j++)
      {
        o_vec[i] += this->GetIJ(i, j) * i_vec[j];
      }
    }
    Vector3 oV(o_vec[0], o_vec[1], o_vec[2]);
    return oV;
  }

  /**Set value at row i and column j.*/
  void SetIJ(int i, int j, double value)
  {
    int k = j + 3 * i;
    vals[k] = value;
  }

  /**Add value to value at row i and column j.*/
  void AddIJ(int i, int j, double value)
  {
    int k = j + 3 * i;
    vals[k] += value;
  }

  /**Obtain a copy of the value at row i and column j.*/
  double GetIJ(int i, int j) const
  {
    int k = j + 3 * i;
    return vals[k];
  }

  /**Set row i using a vector.*/
  void SetRowIVec(int i, Vector3 vec)
  {
    vals[0 + 3 * i] = vec.x;
    vals[1 + 3 * i] = vec.y;
    vals[2 + 3 * i] = vec.z;
  }

  /**Set column j using a vector.*/
  void SetColJVec(int j, Vector3 vec)
  {
    vals[j + 3 * 0] = vec.x;
    vals[j + 3 * 1] = vec.y;
    vals[j + 3 * 2] = vec.z;
  }

  /**Sets the diagonal of the matrix.*/
  void SetDiagonalVec(double a00, double a11, double a22)
  {
    vals[0 + 3 * 0] = a00;
    vals[1 + 3 * 1] = a11;
    vals[2 + 3 * 2] = a22;
  }

  /**Get the determinant using specified row [default:0].*/
  double Det(int row = 0)
  {
    double det = 0.0;

    int sign = -1;
    for (int j = 0; j < 3; j++)
    {
      int k = j + 3 * row;
      sign *= -1;

      det += sign * GetIJ(row, j) * MinorIJ(row, j);
    }

    return det;
  }

  /**Get the minor value associated with row ir and column jr.*/
  double MinorIJ(int ir, int jr)
  {
    double a[] = {0.0, 0.0, 0.0, 0.0};

    int k = -1;
    for (int i = 0; i < 3; i++)
    {
      if (i == ir)
        continue;

      for (int j = 0; j < 3; j++)
      {
        if (j == jr)
          continue;

        k++;
        int kr = j + 3 * i;
        a[k] = vals[kr];
      }
    }

    return a[0] * a[3] - a[2] * a[1];
  }

  /**Compute the matrix transpose.*/
  Matrix3x3 Transpose()
  {
    Matrix3x3 oM;
    for (int i = 0; i < 3; i++)
    {
      for (int j = 0; j < 3; j++)
      {
        int kO = j + 3 * i;
        int kI = i + 3 * j;

        oM.vals[kO] = vals[kI];
      }
    }
    return oM;
  }

  /**Compute the matrix transpose.*/
  Matrix3x3 Inverse()
  {
    Matrix3x3 oM;
    Matrix3x3 oMT;

    // Compute matrix of minors
    for (int i = 0; i < 3; i++)
    {
      for (int j = 0; j < 3; j++)
      {
        oM.SetIJ(i, j, MinorIJ(i, j));
      }
    }

    // Compute matrix of cofactors
    int sign = -1;
    for (int i = 0; i < 3; i++)
    {
      for (int j = 0; j < 3; j++)
      {
        sign *= -1;
        oM.SetIJ(i, j, oM.GetIJ(i, j) * sign);
      }
    }

    // Compute the transpose
    oMT = oM.Transpose();

    // Get determinant
    double det = Det();

    return oMT * (1.0 / det);
  }

  /**Outputs the matrix to a stream.*/
  friend std::ostream& operator<<(std::ostream& out, Matrix3x3& inM)
  {
    out << "[";
    for (int i = 0; i < 3; i++)
    {
      for (int j = 0; j < 3; j++)
      {
        out << inM.GetIJ(i, j) << " ";
      }
      if (i != 2)
        out << "\n";
      else
        out << "]\n";
    }

    return out;
  }

  std::string PrintS()
  {
    std::stringstream out;
    out << "[";
    for (int i = 0; i < 3; i++)
    {
      for (int j = 0; j < 3; j++)
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
