// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "framework/math/vector3.h"

namespace opensn
{

// Scalar multiplication for Vector3
Vector3
operator*(const double value, const Vector3& vec)
{
  return Vector3(vec.x * value, vec.y * value, vec.z * value);
}

} // namespace opensn
