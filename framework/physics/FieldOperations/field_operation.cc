#include "framework/physics/FieldOperations/field_operation.h"

namespace chi_physics::field_operations
{

// Since there are no input parameters we will not register this object

// ##################################################################
/**Returns the input parameters.*/
chi::InputParameters
FieldOperation::GetInputParameters()
{
  return ChiObject::GetInputParameters();
}

// ##################################################################
/**Constructor.*/
FieldOperation::FieldOperation(const chi::InputParameters& params) : ChiObject(params)
{
}

} // namespace chi_physics::field_operations
