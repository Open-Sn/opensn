// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "lua/framework/console/console.h"
#include "modules/point_reactor_kinetics/point_reactor_kinetics.h"
#include "framework/runtime.h"
#include "framework/logging/log.h"

using namespace opensn;

namespace opensnlua
{

namespace
{

InputParameters
GetSyntax_SetParam()
{
  InputParameters params;

  params.SetGeneralDescription(
    "Lua wrapper function for setting parameters in the PointReactorKinetics module.");
  params.SetDocGroup("prk");

  params.AddRequiredParameter<size_t>("arg0", "Handle to a <TT>PRKSolver</TT> object.");
  params.AddRequiredParameter<std::string>("arg1", "Text name of the parameter to set.");

  params.AddRequiredParameter<double>("arg2", "Value to set to the parameter pointed to by arg1");

  params.ConstrainParameterRange("arg1", AllowableRangeList::New({"rho"}));

  return params;
}

ParameterBlock
SetParam(const InputParameters& params)
{
  const std::string fname = __FUNCTION__;
  const size_t handle = params.GetParamValue<size_t>("arg0");

  auto& solver = opensn::GetStackItem<opensn::PRKSolver>(opensn::object_stack, handle, fname);

  const auto param_name = params.GetParamValue<std::string>("arg1");
  const auto& value_param = params.GetParam("arg2");

  if (param_name == "rho")
  {
    OpenSnInvalidArgumentIf(value_param.Type() != ParameterBlockType::FLOAT,
                            "If arg1 is \"rho\" then arg2 must be of type FLOAT");
    solver.SetRho(value_param.GetValue<double>());
  }
  else
    OpenSnInvalidArgument("Invalid property name \"" + param_name);

  return ParameterBlock(); // Return empty param block
}

InputParameters
GetParamSyntax()
{
  InputParameters params;

  params.SetGeneralDescription(
    "Lua wrapper function for getting parameters from the PointReactorKinetics"
    " module.");
  params.SetDocGroup("prk");

  params.AddRequiredParameter<size_t>("arg0", "Handle to a <TT>PRKSolver</TT> object.");
  params.AddRequiredParameter<std::string>("arg1", "Text name of the parameter to get.");

  params.ConstrainParameterRange(
    "arg1",
    AllowableRangeList::New(
      {"population_prev", "population_next", "period", "time_prev", "time_next"}));
  return params;
}

ParameterBlock
GetParam(const InputParameters& params)
{
  const std::string fname = __FUNCTION__;
  const size_t handle = params.GetParamValue<size_t>("arg0");

  auto& solver = opensn::GetStackItem<opensn::PRKSolver>(opensn::object_stack, handle, fname);

  const auto param_name = params.GetParamValue<std::string>("arg1");
  ParameterBlock outputs;

  if (param_name == "population_prev")
    outputs.AddParameter("", solver.PopulationPrev());
  else if (param_name == "population_next")
    outputs.AddParameter("", solver.PopulationNew());
  else if (param_name == "period")
    outputs.AddParameter("", solver.Period());
  else if (param_name == "time_prev")
    outputs.AddParameter("", solver.TimePrev());
  else if (param_name == "time_next")
    outputs.AddParameter("", solver.TimeNew());
  else
    OpenSnInvalidArgument("Invalid property name \"" + param_name);

  return outputs;
}

} // namespace

RegisterWrapperFunctionInNamespace(prk, GetParam, GetParamSyntax, GetParam);
RegisterWrapperFunctionInNamespace(prk, SetParam, GetSyntax_SetParam, SetParam);

} // namespace opensnlua
