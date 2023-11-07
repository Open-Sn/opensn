#include "framework/mesh/logical_volume/logical_volume_interface.h"

#include "framework/mesh/logical_volume/logical_volume.h"

#include "framework/runtime.h"

namespace chi_mesh
{

chi::InputParameters
LogicalVolumeInterface::GetInputParameters()
{
  chi::InputParameters params;

  params.AddOptionalParameter("logical_volume", 0, "Handle to a logical_volume.");

  return params;
}

LogicalVolumeInterface::LogicalVolumeInterface(const chi::InputParameters& params)
  : logical_volume_(
      params.ParametersAtAssignment().Has("logical_volume")
        ? Chi::GetStackItemPtrAsType<const LogicalVolume>(
            Chi::object_stack, params.GetParamValue<size_t>("logical_volume"), __FUNCTION__)
        : nullptr)
{
}

const LogicalVolume*
LogicalVolumeInterface::GetLogicalVolume() const
{
  return logical_volume_ ? &(*logical_volume_) : nullptr;
}

} // namespace chi_mesh
