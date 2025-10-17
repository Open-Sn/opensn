// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include "framework/data_types/vector3.h"
#include <functional>

namespace opensn
{

using ScalarSpatialFunction = std::function<double(const Vector3&)>;
using ScalarSpatialMaterialFunction = std::function<double(int, const Vector3&)>;
using ScalarMaterialFunction = std::function<double(double, int)>;

/// Base class for evaluating spatial material functions given a coordinate.
class VectorSpatialFunction : public std::function<std::vector<double>(const Vector3&, std::size_t)>
{
public:
  VectorSpatialFunction() = default;
  VectorSpatialFunction(const std::function<std::vector<double>(const Vector3&, std::size_t)>& src)
    : std::function<std::vector<double>(const Vector3&, std::size_t)>(src)
  {
  }
};

} // namespace opensn
