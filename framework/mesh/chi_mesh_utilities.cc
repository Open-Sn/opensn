#include "chi_mesh.h"

//###################################################################
/**Tensor product of two vectors.
 * \f$ \vec{\vec{T}} = \vec{x} \otimes \vec{y} \f$*/
chi_mesh::TensorRank2Dim3
chi_mesh::Vector3::OTimes(const Vector3& that) const
{
  TensorRank2Dim3 new_t;
  for (int i = 0; i < 3; ++i)
    for (int j = 0; j < 3; ++j)
      new_t[i](j) = this->operator[](i) * that[j];

  return new_t;
}

//###################################################################
/**Dot product of vector and a rank-2 tensor.
 * \f$ \vec{w} = \vec{x} \bullet \vec{\vec{T}} \f$*/
chi_mesh::Vector3
chi_mesh::Vector3::Dot(const chi_mesh::TensorRank2Dim3& that) const
{
  chi_mesh::Vector3 new_vec;
  for (int i = 0; i < 3; ++i)
    new_vec(i) = this->Dot(that.t[i]);

  return new_vec;
}

//###################################################################
/**Returns a 3D vector multiplied by the given scalar from the left.
 * \f$ \vec{w} = \alpha \vec{x}\f$*/
chi_mesh::Vector3
operator*(const double value, const chi_mesh::Vector3& that)
{
  chi_mesh::Vector3 newVector;
  newVector.x = that.x * value;
  newVector.y = that.y * value;
  newVector.z = that.z * value;

  return newVector;
}

//###################################################################
/**Dot product of rank-2 tensor with a vector.
 * \f$ \vec{w} = \vec{\vec{T}} \bullet \vec{x} \f$*/
chi_mesh::Vector3
chi_mesh::TensorRank2Dim3::Dot(const chi_mesh::Vector3& v) const
{
  chi_mesh::Vector3 newVector;
  for (int i = 0; i < 3; ++i)
    newVector(i) = t[i].Dot(v);

  return newVector;
}

//###################################################################
/**Returns the diagonal of a rank-2 dim-3 tensor as a vector3.
 * \f$ \vec{w} = \text{diag} \vec{\vec{T}} \f$*/
chi_mesh::Vector3
chi_mesh::TensorRank2Dim3::Diag() const
{
  chi_mesh::Vector3 newVector;
  for (int i = 0; i < 3; ++i)
    newVector(i) = t[i][i];

  return newVector;
}

//###################################################################
/**Rank-2 dim-3 tensor multiplied from the left with a scalar.
 * \f$ \vec{\vec{W}} = \alpha \vec{\vec{T}} \f$*/
chi_mesh::TensorRank2Dim3
operator*(const double value, const chi_mesh::TensorRank2Dim3& that)
{
  chi_mesh::TensorRank2Dim3 new_t = that;
  for (int i = 0; i < 3; ++i)
  {
    for (int j = 0; j < 3; ++j)
      new_t[i](j) *= value;
  }

  return new_t;
}