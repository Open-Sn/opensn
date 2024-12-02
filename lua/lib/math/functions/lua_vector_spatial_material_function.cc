// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "lua/lib/math/functions/lua_vector_spatial_material_function.h"
#include "framework/runtime.h"
#include "lua/lib/console.h"
#include "lua/lib/types.h"
#include "LuaBridge/LuaBridge.h"
#include <lua.h>

using namespace opensn;

namespace opensnlua
{

OpenSnRegisterObjectInNamespace(opensn, LuaVectorSpatialMaterialFunction);

InputParameters
LuaVectorSpatialMaterialFunction::GetInputParameters()
{
  InputParameters params = VectorSpatialMaterialFunction::GetInputParameters();
  params.AddRequiredParameter<std::string>("function_name", "Name of the lua function");
  return params;
}

std::shared_ptr<LuaVectorSpatialMaterialFunction>
LuaVectorSpatialMaterialFunction::Create(const opensn::ParameterBlock& params)
{
  auto& factory = opensn::ObjectFactory::GetInstance();
  return factory.Create<LuaVectorSpatialMaterialFunction>(
    "opensn::LuaVectorSpatialMaterialFunction", params);
}

LuaVectorSpatialMaterialFunction::LuaVectorSpatialMaterialFunction(const InputParameters& params)
  : opensn::VectorSpatialMaterialFunction(params),
    function_name_(params.GetParamValue<std::string>("function_name"))
{
}

std::vector<double>
LuaVectorSpatialMaterialFunction::Evaluate(const opensn::Vector3& xyz,
                                           int mat_id,
                                           int num_components) const
{
  // Load lua function
  lua_State* L = console.GetConsoleState();
  auto fn = luabridge::getGlobal(L, function_name_.c_str());
  auto lua_return = fn(xyz, mat_id);

  // Check return value
  OpenSnLogicalErrorIf(lua_return.size() != num_components,
                       "Call to lua function " + function_name_ + " returned a vector of size " +
                         std::to_string(lua_return.size()) +
                         ", which is not the same as the number of groups " +
                         std::to_string(num_components) + ".");
  std::vector<double> ret(lua_return.size());
  for (auto i = 0; i < lua_return.size(); i++)
    ret[i] = lua_return[i];
  return ret;
}

} // namespace opensnlua
