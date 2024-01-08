#include "framework/math/functions/scalar_spatial_material_function.h"

namespace opensn
{

InputParameters
ScalarSpatialMaterialFunction::GetInputParameters()
{
  InputParameters params = Function::GetInputParameters();
  return params;
}

ScalarSpatialMaterialFunction::ScalarSpatialMaterialFunction(const InputParameters& params)
  : Function(params)
{
}

} // namespace opensn
