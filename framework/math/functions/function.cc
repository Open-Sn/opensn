#include "function.h"

namespace opensn
{

InputParameters
Function::GetInputParameters()
{
  InputParameters params = Object::GetInputParameters();
  return params;
}

Function::Function(const InputParameters& params) : Object(params)
{
}

} // namespace opensn
