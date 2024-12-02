// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "lua/lib/math/functions/lua_scalar_spatial_material_function.h"
#include "framework/runtime.h"
#include "lua/lib/console.h"
#include "lua/lib/types.h"
#include "LuaBridge/LuaBridge.h"
#include <lua.h>

using namespace opensn;

namespace opensnlua
{

OpenSnRegisterObjectInNamespace(opensn, LuaScalarSpatialMaterialFunction);

InputParameters
LuaScalarSpatialMaterialFunction::GetInputParameters()
{
  InputParameters params = ScalarSpatialMaterialFunction::GetInputParameters();
  params.AddRequiredParameter<std::string>("function_name", "Name of the lua function");
  return params;
}

std::shared_ptr<LuaScalarSpatialMaterialFunction>
LuaScalarSpatialMaterialFunction::Create(const opensn::ParameterBlock& params)
{
  auto& factory = opensn::ObjectFactory::GetInstance();
  return factory.Create<LuaScalarSpatialMaterialFunction>(
    "opensn::LuaScalarSpatialMaterialFunction", params);
}

LuaScalarSpatialMaterialFunction::LuaScalarSpatialMaterialFunction(const InputParameters& params)
  : ScalarSpatialMaterialFunction(params),
    function_name_(params.GetParamValue<std::string>("function_name"))
{
}

double
LuaScalarSpatialMaterialFunction::Evaluate(int mat_id, const opensn::Vector3& xyz) const
{
  lua_State* L = console.GetConsoleState();
  auto fn = luabridge::getGlobal(L, function_name_.c_str());
  auto res = fn(mat_id, xyz);
  return res[0];
}

} // namespace opensnlua
