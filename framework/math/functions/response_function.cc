#include "framework/math/functions/response_function.h"

namespace opensn
{

InputParameters
ResponseFunction::GetInputParameters()
{
  InputParameters params = Function::GetInputParameters();
  return params;
}

ResponseFunction::ResponseFunction(const InputParameters& params) : Function(params)
{
}

} // namespace opensn
