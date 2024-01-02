#include "framework/physics/field_operations/field_operation.h"

namespace opensn
{

// Since there are no input parameters we will not register this object

InputParameters
FieldOperation::GetInputParameters()
{
  return Object::GetInputParameters();
}

FieldOperation::FieldOperation(const InputParameters& params) : Object(params)
{
}

} // namespace opensn
