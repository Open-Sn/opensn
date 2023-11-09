#include "framework/mesh/mesh.h"

namespace chi_mesh
{

Vector3
TensorRank2Dim3::Dot(const Vector3& v) const
{
  Vector3 newVector;
  for (int i = 0; i < 3; ++i)
    newVector(i) = t[i].Dot(v);

  return newVector;
}

Vector3
TensorRank2Dim3::Diag() const
{
  Vector3 newVector;
  for (int i = 0; i < 3; ++i)
    newVector(i) = t[i][i];

  return newVector;
}

TensorRank2Dim3
operator*(const double value, const TensorRank2Dim3& that)
{
  TensorRank2Dim3 new_t = that;
  for (int i = 0; i < 3; ++i)
  {
    for (int j = 0; j < 3; ++j)
      new_t[i](j) *= value;
  }

  return new_t;
}

} // namespace chi_mesh
