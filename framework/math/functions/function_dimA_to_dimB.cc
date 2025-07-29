// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "framework/math/functions/function_dimA_to_dimB.h"
#include "framework/logging/log.h"

namespace opensn
{

InputParameters
FunctionDimAToDimB::GetInputParameters()
{
  InputParameters params;

  params.AddRequiredParameter<size_t>(
    "input_dimension", "The dimension of the input values (excluding the position).");

  params.AddRequiredParameter<size_t>(
    "output_dimension", "The dimension of the output values (excluding the position).");

  return params;
}

FunctionDimAToDimB::FunctionDimAToDimB(const InputParameters& params)
  : input_dimension_(params.GetParamValue<size_t>("input_dimension")),
    output_dimension_(params.GetParamValue<size_t>("output_dimension"))
{
}

double
FunctionDimAToDimB::GetScalarFunction1Parameter(double) const
{
  OpenSnLogicalError("No available function");
}
double
FunctionDimAToDimB::GetScalarFunctionSlope1Parameter(double) const
{
  OpenSnLogicalError("No available function");
}
double
FunctionDimAToDimB::GetScalarFunctionCurvature1Parameter(double) const
{
  OpenSnLogicalError("No available function");
}

double
FunctionDimAToDimB::GetScalarFunction4Parameters(double, double, double, double) const
{
  OpenSnLogicalError("No available function");
}

double
FunctionDimAToDimB::GetScalarFunctionSlope4Parameters(double, double, double, double) const
{
  OpenSnLogicalError("No available function");
}
double
FunctionDimAToDimB::GetScalarFunctionCurvature4Parameters(double, double, double, double) const
{
  OpenSnLogicalError("No available function");
}
} // namespace opensn
