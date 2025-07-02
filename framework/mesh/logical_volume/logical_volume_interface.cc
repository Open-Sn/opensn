// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "framework/mesh/logical_volume/logical_volume_interface.h"
#include "framework/mesh/logical_volume/logical_volume.h"
#include "framework/runtime.h"

namespace opensn
{

InputParameters
LogicalVolumeInterface::GetInputParameters()
{
  InputParameters params;

  params.AddOptionalParameter(
    "logical_volume", std::shared_ptr<LogicalVolume>{}, "Handle to a logical_volume.");

  return params;
}

LogicalVolumeInterface::LogicalVolumeInterface(const InputParameters& params)
  : logical_volume_(params.GetSharedPtrParam<LogicalVolume>("logical_volume", false))
{
}

const std::shared_ptr<LogicalVolume>
LogicalVolumeInterface::GetLogicalVolume() const
{
  return logical_volume_;
}

} // namespace opensn
