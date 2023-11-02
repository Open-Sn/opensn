#include "LogicalVolume.h"

namespace chi_mesh
{

chi::InputParameters LogicalVolume::GetInputParameters()
{
  return ChiObject::GetInputParameters();
}

LogicalVolume::LogicalVolume(const chi::InputParameters& params)
  : ChiObject(params)
{
}

} // namespace chi_mesh