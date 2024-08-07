// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "framework/math/functions/function_dimA_to_dimB.h"
#include "framework/logging/log.h"

namespace opensn
{

InputParameters
FunctionDimAToDimB::GetInputParameters()
{
  InputParameters params = Object::GetInputParameters();

  params.AddRequiredParameter<size_t>(
    "input_dimension", "The dimension of the input values (excluding the position).");

  params.AddRequiredParameter<size_t>(
    "output_dimension", "The dimension of the output values (excluding the position).");

  return params;
}

FunctionDimAToDimB::FunctionDimAToDimB(const InputParameters& params)
  : Object(params),
    input_dimension_(params.GetParamValue<size_t>("input_dimension")),
    output_dimension_(params.GetParamValue<size_t>("output_dimension"))
{
}

double
FunctionDimAToDimB::ScalarFunction1Parameter(double) const
{
  OpenSnLogicalError("No available function");
}
double
FunctionDimAToDimB::ScalarFunctionSlope1Parameter(double) const
{
  OpenSnLogicalError("No available function");
}
double
FunctionDimAToDimB::ScalarFunctionCurvature1Parameter(double) const
{
  OpenSnLogicalError("No available function");
}

double
FunctionDimAToDimB::ScalarFunction4Parameters(double, double, double, double) const
{
  OpenSnLogicalError("No available function");
}

double
FunctionDimAToDimB::ScalarFunctionSlope4Parameters(double, double, double, double) const
{
  OpenSnLogicalError("No available function");
}
double
FunctionDimAToDimB::ScalarFunctionCurvature4Parameters(double, double, double, double) const
{
  OpenSnLogicalError("No available function");
}
} // namespace opensn
