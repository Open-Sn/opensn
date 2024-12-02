// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "lua/lib/math/functions/lua_vector_spatial_function.h"
#include "framework/runtime.h"
#include "lua/lib/console.h"
#include "lua/lib/types.h"
#include "LuaBridge/LuaBridge.h"
#include <lua.h>

using namespace opensn;

namespace opensnlua
{

OpenSnRegisterObjectInNamespace(opensn, LuaVectorSpatialFunction);

InputParameters
LuaVectorSpatialFunction::GetInputParameters()
{
  InputParameters params = VectorSpatialFunction::GetInputParameters();
  params.AddRequiredParameter<std::string>("function_name", "Name of the lua function");
  return params;
}

std::shared_ptr<LuaVectorSpatialFunction>
LuaVectorSpatialFunction::Create(const opensn::ParameterBlock& params)
{
  auto& factory = opensn::ObjectFactory::GetInstance();
  return factory.Create<LuaVectorSpatialFunction>("opensn::LuaVectorSpatialFunction", params);
}

LuaVectorSpatialFunction::LuaVectorSpatialFunction(const opensn::InputParameters& params)
  : opensn::VectorSpatialFunction(params),
    function_name_(params.GetParamValue<std::string>("function_name"))
{
}

std::vector<double>
LuaVectorSpatialFunction::Evaluate(const opensn::Vector3& xyz, int num_components) const
{
  // Load lua function
  lua_State* L = console.GetConsoleState();
  auto fn = luabridge::getGlobal(L, function_name_.c_str());
  auto lua_return = fn(xyz);
  OpenSnLogicalErrorIf(lua_return.hasFailed(), "Call to " + function_name_ + " failed.");
  auto tab = lua_return[0];
  OpenSnLogicalErrorIf(not tab.isTable(),
                       "The return value from " + function_name_ + " must be a table.");
  // Check return value
  OpenSnLogicalErrorIf(tab.length() != num_components,
                       "Call to lua function " + function_name_ + " returned a vector of size " +
                         std::to_string(tab.length()) +
                         ", which is not the same as the number of groups " +
                         std::to_string(num_components) + ".");
  std::vector<double> ret(tab.length());
  for (auto i = 0; i < tab.length(); i++)
    ret[i] = tab[i + 1];
  return ret;
}

} // namespace opensnlua
