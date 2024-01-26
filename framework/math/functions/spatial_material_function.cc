#include "framework/math/functions/spatial_material_function.h"

namespace opensn
{

InputParameters
SpatialMaterialFunction::GetInputParameters()
{
  InputParameters params = Function::GetInputParameters();
  return params;
}

SpatialMaterialFunction::SpatialMaterialFunction(const InputParameters& params) : Function(params)
{
}

} // namespace opensn
