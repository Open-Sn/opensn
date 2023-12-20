#include "framework/mesh/mesh.h"

namespace opensn
{

Vector3
operator*(const double value, const Vector3& that)
{
  Vector3 newVector;
  newVector.x = that.x * value;
  newVector.y = that.y * value;
  newVector.z = that.z * value;

  return newVector;
}

} // namespace opensn
