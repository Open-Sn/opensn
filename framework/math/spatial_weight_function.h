// SPDX-FileCopyrightText: 2025 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include "framework/data_types/vector3.h"
#include "framework/mesh/mesh/mesh.h"

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

/**
 * Factor that converts an integral computed with the weights above into a physical integral.
 *
 * The cylindrical and spherical weights omit the angular extent of the coordinate system, so
 * integrals built from them are per radian in RZ and per steradian in 1D spherical geometry.
 * The physical measure is 2*pi*r in RZ and 4*pi*r^2 in 1D spherical geometry.
 *
 * \return 2*pi for cylindrical, 4*pi for spherical, and 1 for Cartesian coordinates.
 */
double IntegralMeasureScale(CoordinateSystemType coord_sys);

} // namespace opensn
