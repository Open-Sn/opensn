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

  params.AddOptionalParameter("logical_volume", 0, "Handle to a logical_volume.");

  return params;
}

LogicalVolumeInterface::LogicalVolumeInterface(const InputParameters& params)
  : logical_volume_(
      params.ParametersAtAssignment().Has("logical_volume")
        ? GetStackItemPtrAsType<const LogicalVolume>(
            object_stack, params.GetParamValue<size_t>("logical_volume"), __FUNCTION__)
        : nullptr)
{
}

const LogicalVolume*
LogicalVolumeInterface::GetLogicalVolume() const
{
  return logical_volume_ ? &(*logical_volume_) : nullptr;
}

} // namespace opensn
