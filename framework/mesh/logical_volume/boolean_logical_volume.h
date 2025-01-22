// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include "framework/mesh/logical_volume/logical_volume.h"

namespace opensn
{

/// Boolean volume
class BooleanLogicalVolume : public LogicalVolume
{
public:
  std::vector<std::pair<bool, std::shared_ptr<const LogicalVolume>>> parts;

  explicit BooleanLogicalVolume(const InputParameters& params);

  bool Inside(const Vector3& point) const override;

public:
  static InputParameters GetInputParameters();
  static std::shared_ptr<BooleanLogicalVolume> Create(const ParameterBlock& params);
};

} // namespace opensn
