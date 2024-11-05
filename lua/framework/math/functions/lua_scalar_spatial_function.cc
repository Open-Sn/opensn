// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "lua/framework/math/functions/lua_scalar_spatial_function.h"
#include "lua/framework/lua.h"
#include "lua/framework/console/console.h"
#include "framework/runtime.h"
#include "framework/object_factory.h"

using namespace opensn;

namespace opensnlua
{

OpenSnRegisterObjectInNamespace(opensn, LuaScalarSpatialFunction);

InputParameters
LuaScalarSpatialFunction::GetInputParameters()
{
  InputParameters params = ScalarSpatialFunction::GetInputParameters();
  params.AddRequiredParameter<std::string>("lua_function_name", "Name of the lua function");
  return params;
}

LuaScalarSpatialFunction::LuaScalarSpatialFunction(const InputParameters& params)
  : ScalarSpatialFunction(params),
    lua_function_name_(params.ParamValue<std::string>("lua_function_name"))
{
}

double
LuaScalarSpatialFunction::Evaluate(const opensn::Vector3& xyz) const
{
  lua_State* L = console.GetConsoleState();
  return LuaCall<double>(L, lua_function_name_, xyz);
}

} // namespace opensnlua
