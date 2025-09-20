// SPDX-FileCopyrightText: 2025 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include "framework/data_types/vector3.h"
#include "framework/mesh/mesh.h"

namespace opensn
{

struct SpatialWeightFunction
{
  virtual double operator()(const Vector3& pt) const = 0;
  virtual ~SpatialWeightFunction() = default;

  static std::shared_ptr<SpatialWeightFunction> FromCoordinateType(CoordinateSystemType coord_sys);
};

struct CartesianSpatialWeightFunction : public SpatialWeightFunction
{
  double operator()(const Vector3& pt) const override { return 1.0; }
};

struct SphericalSpatialWeightFunction : public CartesianSpatialWeightFunction
{
  double operator()(const Vector3& pt) const override { return pt[2] * pt[2]; }
};

struct CylindricalSpatialWeightFunction : public CartesianSpatialWeightFunction
{
  double operator()(const Vector3& pt) const override { return pt[0]; }
};

} // namespace opensn
