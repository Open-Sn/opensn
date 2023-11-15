#include "framework/math/functions/scalar_material_function.h"

namespace opensn
{

InputParameters
ScalarMaterialFunction::GetInputParameters()
{
  InputParameters params = Function::GetInputParameters();
  return params;
}

ScalarMaterialFunction::ScalarMaterialFunction(const InputParameters& params) : Function(params)
{
}

} // namespace opensn
