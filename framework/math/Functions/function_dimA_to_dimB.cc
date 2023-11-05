#include "opensn/framework/math/Functions/function_dimA_to_dimB.h"

#include "opensn/framework/logging/chi_log.h"

namespace chi_math
{

chi::InputParameters
FunctionDimAToDimB::GetInputParameters()
{
  chi::InputParameters params = ChiObject::GetInputParameters();

  params.AddRequiredParameter<size_t>(
    "input_dimension", "The dimension of the input values (excluding the position).");

  params.AddRequiredParameter<size_t>(
    "output_dimension", "The dimension of the output values (excluding the position).");

  return params;
}

FunctionDimAToDimB::FunctionDimAToDimB(const chi::InputParameters& params)
  : ChiObject(params),
    input_dimension_(params.GetParamValue<size_t>("input_dimension")),
    output_dimension_(params.GetParamValue<size_t>("output_dimension"))
{
}

double
FunctionDimAToDimB::ScalarFunction1Parameter(double) const
{
  ChiLogicalError("No available function");
}
double
FunctionDimAToDimB::ScalarFunctionSlope1Parameter(double) const
{
  ChiLogicalError("No available function");
}
double
FunctionDimAToDimB::ScalarFunctionCurvature1Parameter(double) const
{
  ChiLogicalError("No available function");
}

double
FunctionDimAToDimB::ScalarFunction4Parameters(double, double, double, double) const
{
  ChiLogicalError("No available function");
}

double
FunctionDimAToDimB::ScalarFunctionSlope4Parameters(double, double, double, double) const
{
  ChiLogicalError("No available function");
}
double
FunctionDimAToDimB::ScalarFunctionCurvature4Parameters(double, double, double, double) const
{
  ChiLogicalError("No available function");
}
} // namespace chi_math
