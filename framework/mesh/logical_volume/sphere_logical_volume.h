// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include "framework/mesh/logical_volume/logical_volume.h"

namespace opensn
{

/**Spherical logical volume.*/
class SphereLogicalVolume : public LogicalVolume
{
public:
  static InputParameters GetInputParameters();
  explicit SphereLogicalVolume(const InputParameters& params);

  bool Inside(const Vector3& point) const override;

protected:
  double r_;
  double x0_, y0_, z0_;
};

} // namespace opensn
