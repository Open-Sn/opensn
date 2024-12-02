// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "lua/lib/math/functions/lua_scalar_spatial_function.h"
#include "framework/runtime.h"
#include "lua/lib/console.h"
#include "lua/lib/types.h"
#include "LuaBridge/LuaBridge.h"
#include <lua.h>

using namespace opensn;

namespace opensnlua
{

OpenSnRegisterObjectInNamespace(opensn, LuaScalarSpatialFunction);

InputParameters
LuaScalarSpatialFunction::GetInputParameters()
{
  InputParameters params = ScalarSpatialFunction::GetInputParameters();
  params.AddRequiredParameter<std::string>("function_name", "Name of the lua function");
  return params;
}

std::shared_ptr<LuaScalarSpatialFunction>
LuaScalarSpatialFunction::Create(const opensn::ParameterBlock& params)
{
  auto& factory = opensn::ObjectFactory::GetInstance();
  return factory.Create<LuaScalarSpatialFunction>("opensn::LuaScalarSpatialFunction", params);
}

LuaScalarSpatialFunction::LuaScalarSpatialFunction(const InputParameters& params)
  : ScalarSpatialFunction(params),
    function_name_(params.GetParamValue<std::string>("function_name"))
{
}

double
LuaScalarSpatialFunction::Evaluate(const opensn::Vector3& xyz) const
{
  lua_State* L = console.GetConsoleState();
  auto fn = luabridge::getGlobal(L, function_name_.c_str());
  return fn(xyz)[0];
}

} // namespace opensnlua
