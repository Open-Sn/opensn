// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include "framework/parameters/input_parameters.h"

namespace opensn
{

class LogicalVolume;

/**Interface class to add a dependency on a logical volume. Two things need to
 * be done to use this interface. 1) Derive from it. 2) Add its parameters to
 * the child class. Now it will require a handle to logical volume in the input
 * language.*/
class LogicalVolumeInterface
{
protected:
  static InputParameters GetInputParameters();

  explicit LogicalVolumeInterface(const InputParameters& params);

  const LogicalVolume* GetLogicalVolume() const;

private:
  const std::shared_ptr<const LogicalVolume> logical_volume_;
};

} // namespace opensn
