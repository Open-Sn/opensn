#include "framework/lua.h"
#include "framework/console/console.h"

#include "modules/point_reactor_kinetics/point_reactor_kinetics.h"
#include "prk_lua_utils.h"

#include "framework/runtime.h"
#include "framework/logging/log.h"

using namespace opensn;

namespace opensnlua::prk
{

RegisterLuaFunctionNamespace(PRKGetParam, prk, GetParam);
RegisterLuaFunctionNamespace(PRKSetParam, prk, SetParam);

int
PRKGetParam(lua_State* L)
{
  const std::string fname = "prk.GetParam";
  LuaCheckArgs<int, std::string>(L, fname);

  const auto handle = LuaArg<int>(L, 1);
  auto solver =
    opensn::GetStackItem<opensn::prk::TransientSolver>(opensn::object_stack, handle, fname);

  const auto param_name = LuaArg<std::string>(L, 2);
  if (param_name == "population_prev")
    return LuaReturn(L, solver.PopulationPrev());
  else if (param_name == "population_next")
    return LuaReturn(L, solver.PopulationNew());
  else if (param_name == "period")
    return LuaReturn(L, solver.Period());
  else if (param_name == "time_prev")
    return LuaReturn(L, solver.TimePrev());
  else if (param_name == "time_next")
    return LuaReturn(L, solver.TimeNew());
  else
    throw std::invalid_argument(fname + ": Invalid parameter '" + param_name + "'.");
}

int
PRKSetParam(lua_State* L)
{
  const std::string fname = "prk.SetParam";
  LuaCheckArgs<size_t, std::string, double>(L, fname);

  const auto handle = LuaArg<size_t>(L, 1);
  auto& solver =
    opensn::GetStackItem<opensn::prk::TransientSolver>(opensn::object_stack, handle, fname);

  const auto param_name = LuaArg<std::string>(L, 2);
  if (param_name == "rho")
  {
    const auto val = LuaArg<double>(L, 3);
    solver.SetRho(val);
  }
  else
    throw std::invalid_argument(fname + ": Invalid parameter '" + param_name + "'\"'.");

  return LuaReturn(L);
}

RegisterWrapperFunctionNamespace(prk, SetParam, GetSyntax_SetParam, SetParam);

InputParameters
GetSyntax_SetParam()
{
  InputParameters params;

  params.SetGeneralDescription(
    "Lua wrapper function for setting parameters in the PointReactorKinetics"
    " module.");
  params.SetDocGroup("prk");

  params.AddRequiredParameter<size_t>("arg0", "Handle to a <TT>prk::TransientSolver</TT> object.");
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

  auto& solver =
    opensn::GetStackItem<opensn::prk::TransientSolver>(opensn::object_stack, handle, fname);

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

RegisterWrapperFunctionNamespace(prk, GetParam, GetParamSyntax, GetParam);

InputParameters
GetParamSyntax()
{
  InputParameters params;

  params.SetGeneralDescription(
    "Lua wrapper function for getting parameters from the PointReactorKinetics"
    " module.");
  params.SetDocGroup("prk");

  params.AddRequiredParameter<size_t>("arg0", "Handle to a <TT>prk::TransientSolver</TT> object.");
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

  auto& solver =
    opensn::GetStackItem<opensn::prk::TransientSolver>(opensn::object_stack, handle, fname);

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

} // namespace opensnlua::prk
