// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "physics/problems/problem.h"
#include "physics/solvers/time_steppers/time_stepper.h"
#include "framework/runtime.h"
#include "framework/logging/log.h"
#include "physics/solvers/time_steppers/constant_time_stepper.h"
#include "framework/object_factory.h"

namespace opensn
{

InputParameters
Problem::GetInputParameters()
{
  InputParameters params = Object::GetInputParameters();

  params.AddRequiredParameter<std::string>(
    "name",
    "A text name to associate with the solver. This name will be used "
    "in status messages and verbose iterative convergence monitors.");

  return params;
}

Problem::Problem(std::string name) : name_(std::move(name))
{
}

Problem::Problem(const InputParameters& params)
  : Object(params), name_(params.GetParamValue<std::string>("name"))
{
}

std::string
Problem::GetName() const
{
  return name_;
}

std::vector<std::shared_ptr<FieldFunctionGridBased>>&
Problem::GetFieldFunctions()
{
  return field_functions_;
}

const std::vector<std::shared_ptr<FieldFunctionGridBased>>&
Problem::GetFieldFunctions() const
{
  return field_functions_;
}

void
Problem::Initialize()
{
  log.Log() << "\"Initialize()\" method not defined for " << GetName();
}

ParameterBlock
Problem::GetInfo(const ParameterBlock& params) const
{
  return ParameterBlock{};
}

ParameterBlock
Problem::GetInfoWithPreCheck(const ParameterBlock& params) const
{
  if (not params.Has("name"))
  {
    log.LogAllWarning() << "Problem::GetInfo called without "
                           "\"name\" in the parameter list";
    return ParameterBlock{};
  }
  return GetInfo(params);
}

void
Problem::SetProperties(const ParameterBlock& params)
{
}

} // namespace opensn
