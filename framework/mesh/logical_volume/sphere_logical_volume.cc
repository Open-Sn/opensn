// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "framework/mesh/logical_volume/sphere_logical_volume.h"
#include "framework/parameters/input_parameters.h"
#include <cmath>

namespace opensn
{

InputParameters
SphereLogicalVolume::GetInputParameters()
{
  InputParameters params;

  params.AddOptionalParameter("r", 1.0, "Radius of the sphere.");
  params.AddOptionalParameter("x", 0.0, "X-location of the volume.");
  params.AddOptionalParameter("y", 0.0, "Y-location of the volume.");
  params.AddOptionalParameter("z", 0.0, "Z-location of the volume.");

  params.ConstrainParameterRange("r", AllowableRangeLowLimit::New(0.0, false));

  return params;
}

std::shared_ptr<SphereLogicalVolume>
SphereLogicalVolume::Create(const ParameterBlock& params)
{
  return CreateObject<SphereLogicalVolume>("logvol::SphereLogicalVolume", params);
}

SphereLogicalVolume::SphereLogicalVolume(const InputParameters& params)
  : LogicalVolume(params),
    r_(params.GetParamValue<double>("r")),
    x0_(params.GetParamValue<double>("x")),
    y0_(params.GetParamValue<double>("y")),
    z0_(params.GetParamValue<double>("z"))
{
}

bool
SphereLogicalVolume::Inside(const Vector3& point) const
{
  double dx = point.x - x0_;
  double dy = point.y - y0_;
  double dz = point.z - z0_;

  double R2 = dx * dx + dy * dy + dz * dz;

  return R2 <= (r_ * r_);
}

bool
SphereLogicalVolume::GetBoundingBox(Vector3& lo, Vector3& hi) const
{
  // The padding covers the rounding in Inside().
  const double half_width = r_ + 1.0e-12 * (r_ + std::abs(x0_) + std::abs(y0_) + std::abs(z0_));
  lo = Vector3(x0_ - half_width, y0_ - half_width, z0_ - half_width);
  hi = Vector3(x0_ + half_width, y0_ + half_width, z0_ + half_width);
  return true;
}

} // namespace opensn
