#include "framework/mesh/logical_volume/logical_volume.h"

namespace opensn
{

InputParameters
LogicalVolume::GetInputParameters()
{
  return Object::GetInputParameters();
}

LogicalVolume::LogicalVolume(const InputParameters& params) : Object(params)
{
}

} // namespace opensn
