// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "framework/math/time_integrations/time_integration.h"
#include "framework/logging/log_exceptions.h"

namespace opensn
{

InputParameters
TimeIntegration::GetInputParameters()
{
  InputParameters params = Object::GetInputParameters();

  params.AddRequiredParameter<int>("method", "Integer representing time stepping scheme");

  return params;
}

TimeIntegration::TimeIntegration(const InputParameters& params) : Object(params)
{
  const int method_option = params.ParamValue<int>("method");
  if (method_option == static_cast<int>(SteppingMethod::EXPLICIT_EULER))
    method_ = SteppingMethod::EXPLICIT_EULER;
  else if (method_option == static_cast<int>(SteppingMethod::IMPLICIT_EULER))
    method_ = SteppingMethod::IMPLICIT_EULER;
  else if (method_option == static_cast<int>(SteppingMethod::CRANK_NICOLSON))
    method_ = SteppingMethod::CRANK_NICOLSON;
  else if (method_option == static_cast<int>(SteppingMethod::THETA_SCHEME))
    method_ = SteppingMethod::THETA_SCHEME;
  else
    OpenSnInvalidArgument("Unsupported Time Integration scheme");
}

SteppingMethod
TimeIntegration::Method() const
{
  return method_;
}

} // namespace opensn
