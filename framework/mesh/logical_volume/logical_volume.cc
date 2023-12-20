#include "framework/mesh/logical_volume/logical_volume.h"

namespace opensn
{

InputParameters
LogicalVolume::GetInputParameters()
{
  return ChiObject::GetInputParameters();
}

LogicalVolume::LogicalVolume(const InputParameters& params) : ChiObject(params)
{
}

} // namespace opensn
