// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include "framework/data_types/vector3.h"
#include <functional>
#include <vector>

namespace opensn
{

using ScalarSpatialFunction = std::function<double(const Vector3&)>;
using ScalarSpatialMaterialFunction = std::function<double(int, const Vector3&)>;
using ScalarMaterialFunction = std::function<double(double, unsigned int)>;

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

/// Base class for evaluating incoming angular fluxes given group and direction indices.
class AngularFluxFunction : public std::function<double(int, int)>
{
public:
  AngularFluxFunction() = default;
  AngularFluxFunction(const std::function<double(int, int)>& src)
    : std::function<double(int, int)>(src)
  {
  }
};

/// Base class for evaluating group-wise strengths given group index and time.
class GroupTimeFunction : public std::function<double(unsigned int, double)>
{
public:
  GroupTimeFunction() = default;
  GroupTimeFunction(const std::function<double(unsigned int, double)>& src)
    : std::function<double(unsigned int, double)>(src)
  {
  }
};

} // namespace opensn
