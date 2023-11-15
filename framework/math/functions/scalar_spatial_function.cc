#include "framework/math/functions/scalar_spatial_function.h"

namespace opensn
{

InputParameters
ScalarSpatialFunction::GetInputParameters()
{
  InputParameters params = Function::GetInputParameters();
  return params;
}

ScalarSpatialFunction::ScalarSpatialFunction(const InputParameters& params) : Function(params)
{
}

} // namespace opensn
