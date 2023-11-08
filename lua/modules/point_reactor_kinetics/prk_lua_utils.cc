#include "framework/chi_lua.h"
#include "framework/console/chi_console.h"

#include "modules/PointReactorKinetics/point_reactor_kinetics.h"
#include "prk_lua_utils.h"

#include "framework/chi_runtime.h"
#include "framework/logging/chi_log.h"

namespace prk::lua_utils
{

RegisterLuaFunctionAsIs(chiPRKGetParam);
RegisterLuaFunctionAsIs(chiPRKSetParam);

int
chiPRKGetParam(lua_State* L)
{
  const std::string fname = __FUNCTION__;
  const int num_args = lua_gettop(L);
  if (num_args != 2) LuaPostArgAmountError(fname, 2, num_args);

  LuaCheckNilValue(fname, L, 1);
  LuaCheckStringValue(fname, L, 2);

  const int handle = lua_tointeger(L, 1);

  auto solver = Chi::GetStackItem<TransientSolver>(Chi::object_stack, handle, fname);

  const std::string param_name = lua_tostring(L, 2);

  if (param_name == "population_prev") lua_pushnumber(L, solver.PopulationPrev());
  else if (param_name == "population_next")
    lua_pushnumber(L, solver.PopulationNew());
  else if (param_name == "period")
    lua_pushnumber(L, solver.Period());
  else if (param_name == "time_prev")
    lua_pushnumber(L, solver.TimePrev());
  else if (param_name == "time_next")
    lua_pushnumber(L, solver.TimeNew());
  else
    throw std::invalid_argument(fname + ": Invalid parameter \"" + param_name + "\".");

  return 1;
}

int
chiPRKSetParam(lua_State* L)
{
  const std::string fname = __FUNCTION__;
  const int num_args = lua_gettop(L);
  if (num_args != 3) LuaPostArgAmountError(fname, 3, num_args);

  LuaCheckNilValue(fname, L, 1);
  LuaCheckStringValue(fname, L, 2);
  LuaCheckNilValue(fname, L, 3);

  const int handle = lua_tointeger(L, 1);

  auto& solver = Chi::GetStackItem<TransientSolver>(Chi::object_stack, handle, fname);

  const std::string param_name = lua_tostring(L, 2);

  if (param_name == "rho")
  {
    LuaCheckNumberValue(fname + "(handle,\"rho\", : Expects a number value.", L, 3);
    const double val = lua_tonumber(L, 3);
    solver.SetRho(val);
  }
  else
    throw std::invalid_argument(fname + ": Invalid parameter \"" + param_name + "\".");

  return 0;
}

RegisterWrapperFunction(prk, SetParam, GetSyntax_SetParam, SetParam);

chi::InputParameters
GetSyntax_SetParam()
{
  chi::InputParameters params;

  params.SetGeneralDescription(
    "Lua wrapper function for setting parameters in the PointReactorKinetics"
    " module.");
  params.SetDocGroup("prk");

  params.AddRequiredParameter<size_t>("arg0", "Handle to a <TT>prk::TransientSolver</TT> object.");
  params.AddRequiredParameter<std::string>("arg1", "Text name of the parameter to set.");

  params.AddRequiredParameter<double>("arg2", "Value to set to the parameter pointed to by arg1");

  using namespace chi_data_types;
  params.ConstrainParameterRange("arg1", AllowableRangeList::New({"rho"}));

  return params;
}

chi::ParameterBlock
SetParam(const chi::InputParameters& params)
{
  const std::string fname = __FUNCTION__;
  const size_t handle = params.GetParamValue<size_t>("arg0");

  auto& solver = Chi::GetStackItem<TransientSolver>(Chi::object_stack, handle, fname);

  const auto param_name = params.GetParamValue<std::string>("arg1");
  const auto& value_param = params.GetParam("arg2");

  using namespace chi;
  if (param_name == "rho")
  {
    ChiInvalidArgumentIf(value_param.Type() != ParameterBlockType::FLOAT,
                         "If arg1 is \"rho\" then arg2 must be of type FLOAT");
    solver.SetRho(value_param.GetValue<double>());
  }
  else
    ChiInvalidArgument("Invalid property name \"" + param_name);

  return chi::ParameterBlock(); // Return empty param block
}

RegisterWrapperFunction(prk, GetParam, GetParamSyntax, GetParam);

chi::InputParameters
GetParamSyntax()
{
  chi::InputParameters params;

  params.SetGeneralDescription(
    "Lua wrapper function for getting parameters from the PointReactorKinetics"
    " module.");
  params.SetDocGroup("prk");

  params.AddRequiredParameter<size_t>("arg0", "Handle to a <TT>prk::TransientSolver</TT> object.");
  params.AddRequiredParameter<std::string>("arg1", "Text name of the parameter to get.");

  // clang-format off
  using namespace chi_data_types;
  params.ConstrainParameterRange("arg1", AllowableRangeList::New({
    "population_prev", "population_next", "period", "time_prev", "time_next"}));
  // clang-format on
  return params;
}

chi::ParameterBlock
GetParam(const chi::InputParameters& params)
{
  const std::string fname = __FUNCTION__;
  const size_t handle = params.GetParamValue<size_t>("arg0");

  auto& solver = Chi::GetStackItem<TransientSolver>(Chi::object_stack, handle, fname);

  const auto param_name = params.GetParamValue<std::string>("arg1");
  chi::ParameterBlock outputs;

  if (param_name == "population_prev") outputs.AddParameter("", solver.PopulationPrev());
  else if (param_name == "population_next")
    outputs.AddParameter("", solver.PopulationNew());
  else if (param_name == "period")
    outputs.AddParameter("", solver.Period());
  else if (param_name == "time_prev")
    outputs.AddParameter("", solver.TimePrev());
  else if (param_name == "time_next")
    outputs.AddParameter("", solver.TimeNew());
  else
    ChiInvalidArgument("Invalid property name \"" + param_name);

  return outputs;
}

} // namespace prk::lua_utils
