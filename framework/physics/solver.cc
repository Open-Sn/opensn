// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "framework/physics/solver.h"
#include "framework/runtime.h"
#include "framework/logging/log.h"
#include "framework/object_factory.h"

namespace opensn
{

InputParameters
Solver::GetInputParameters()
{
  InputParameters params;

  params.AddRequiredParameter<std::string>(
    "name",
    "A text name to associate with the solver. This name will be used "
    "in status messages and verbose iterative convergence monitors.");

  return params;
}

Solver::Solver(std::string name) : name_(std::move(name))
{
}

Solver::Solver(const InputParameters& params) : name_(params.GetParamValue<std::string>("name"))
{
}

std::string
Solver::GetName() const
{
  return name_;
}

void
Solver::Initialize()
{
  log.Log() << "\"Initialize()\" method not defined for " << GetName();
}

void
Solver::Execute()
{
  log.Log() << "\"Execute()\" method not defined for " << GetName();
}

void
Solver::Step()
{
  log.Log() << "\"Step()\" method not defined for " << GetName();
}

void
Solver::Advance()
{
  log.Log() << "\"Advance()\" method not defined for " << GetName();
}

ParameterBlock
Solver::GetInfo(const ParameterBlock& params) const
{
  return ParameterBlock{};
}

ParameterBlock
Solver::GetInfoWithPreCheck(const ParameterBlock& params) const
{
  if (not params.Has("name"))
  {
    log.LogAllWarning() << "Solver::GetInfo called without "
                           "\"name\" in the parameter list";
    return ParameterBlock{};
  }
  return GetInfo(params);
}

} // namespace opensn
