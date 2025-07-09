// SPDX-FileCopyrightText: 2025 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include "framework/math/vector3.h"
#include "framework/math/geometry.h"

namespace opensn
{

// Define spatial weighting functions
struct SpatialWeightFunction
{
  virtual double operator()(const Vector3& pt) const { return 1.0; }
  virtual ~SpatialWeightFunction() = default;

  static std::shared_ptr<SpatialWeightFunction> FromGeometryType(GeometryType geometry_type);
};

struct SphericalWeightFunction : public SpatialWeightFunction
{
  double operator()(const Vector3& pt) const override { return pt[2] * pt[2]; }
};

struct CylindricalWeightFunction : public SpatialWeightFunction
{
  double operator()(const Vector3& pt) const override { return pt[0]; }
};

} // namespace opensn
