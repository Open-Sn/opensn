// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include "framework/mesh/logical_volume/logical_volume.h"

namespace opensn
{

/**
 * Right Circular Cylinder (RCC) logical volume.
 *
 * Determining whether a point is within an RCC is tricky.
 */
class RCCLogicalVolume : public LogicalVolume
{
public:
  explicit RCCLogicalVolume(const InputParameters& params);

  bool Inside(const Vector3& point) const override;

protected:
  double r_;
  double x0_, y0_, z0_;
  double vx_, vy_, vz_;

public:
  static InputParameters GetInputParameters();
  static std::shared_ptr<RCCLogicalVolume> Create(const ParameterBlock& params);
};

} // namespace opensn
