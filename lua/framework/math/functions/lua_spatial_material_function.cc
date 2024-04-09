// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "lua/framework/math/functions/lua_spatial_material_function.h"
#include "framework/lua.h"
#include "framework/runtime.h"
#include "framework/console/console.h"
#include "framework/object_factory.h"

using namespace opensn;

namespace opensnlua
{

OpenSnRegisterObjectInNamespace(opensn, LuaSpatialMaterialFunction);

InputParameters
LuaSpatialMaterialFunction::GetInputParameters()
{
  InputParameters params = SpatialMaterialFunction::GetInputParameters();
  params.AddRequiredParameter<std::string>("lua_function_name", "Name of the lua function");
  return params;
}

LuaSpatialMaterialFunction::LuaSpatialMaterialFunction(const InputParameters& params)
  : opensn::SpatialMaterialFunction(params),
    lua_function_name_(params.GetParamValue<std::string>("lua_function_name"))
{
}

std::vector<double>
LuaSpatialMaterialFunction::Evaluate(const opensn::Vector3& xyz,
                                     int mat_id,
                                     int num_components) const
{
  // Check response function given
  // Return default if none provided
  if (lua_function_name_.empty())
    return std::vector<double>(num_components, 1.0);

  // Load lua function
  lua_State* L = console.GetConsoleState();
  auto lua_return = LuaCall<std::vector<double>>(L, lua_function_name_, xyz, mat_id);
  // Check return value
  OpenSnLogicalErrorIf(lua_return.size() != num_components,
                       "Call to lua function " + lua_function_name_ +
                         " returned a vector of size " + std::to_string(lua_return.size()) +
                         ", which is not the same as the number of groups " +
                         std::to_string(num_components) + ".");

  return lua_return;
}

} // namespace opensnlua
