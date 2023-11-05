#include "opensn/framework/physics/FieldFunction/GridBasedFieldFunctionInterface.h"

#include "opensn/framework/physics/FieldFunction/fieldfunction_gridbased.h"

namespace chi_physics
{

chi::InputParameters
GridBasedFieldFunctionInterface::GetInputParameters()
{
  chi::InputParameters params = FieldFunctionInterface::GetInputParameters();

  return params;
}

GridBasedFieldFunctionInterface::GridBasedFieldFunctionInterface(const chi::InputParameters& params)
  : FieldFunctionInterface(params)
{
}

FieldFunctionGridBased*
GridBasedFieldFunctionInterface::GetGridBasedFieldFunction() const
{
  auto* ff_ptr = GetFieldFunction();

  auto* grid_based_ff_ptr = dynamic_cast<FieldFunctionGridBased*>(ff_ptr);

  return grid_based_ff_ptr ? grid_based_ff_ptr : nullptr;
}

} // namespace chi_physics
