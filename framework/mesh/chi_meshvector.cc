#include "framework/mesh/chi_mesh.h"

namespace chi_mesh
{

TensorRank2Dim3
Vector3::OTimes(const Vector3& that) const
{
  TensorRank2Dim3 new_t;
  for (int i = 0; i < 3; ++i)
    for (int j = 0; j < 3; ++j)
      new_t[i](j) = this->operator[](i) * that[j];

  return new_t;
}

Vector3
Vector3::Dot(const TensorRank2Dim3& that) const
{
  Vector3 new_vec;
  for (int i = 0; i < 3; ++i)
    new_vec(i) = this->Dot(that.t[i]);

  return new_vec;
}

Vector3
operator*(const double value, const Vector3& that)
{
  Vector3 newVector;
  newVector.x = that.x * value;
  newVector.y = that.y * value;
  newVector.z = that.z * value;

  return newVector;
}

} // namespace chi_mesh
