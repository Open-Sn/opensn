// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "lua/framework/math/functions/lua_scalar_spatial_material_function.h"
#include "framework/lua.h"
#include "framework/runtime.h"
#include "framework/console/console.h"
#include "framework/object_factory.h"

using namespace opensn;

namespace opensnlua
{

OpenSnRegisterObjectInNamespace(opensn, LuaScalarSpatialMaterialFunction);

InputParameters
LuaScalarSpatialMaterialFunction::GetInputParameters()
{
  InputParameters params = ScalarSpatialMaterialFunction::GetInputParameters();
  params.AddRequiredParameter<std::string>("lua_function_name", "Name of the lua function");
  return params;
}

LuaScalarSpatialMaterialFunction::LuaScalarSpatialMaterialFunction(const InputParameters& params)
  : ScalarSpatialMaterialFunction(params),
    lua_function_name_(params.GetParamValue<std::string>("lua_function_name"))
{
}

double
LuaScalarSpatialMaterialFunction::Evaluate(int mat_id, const opensn::Vector3& xyz) const
{
  lua_State* L = console.GetConsoleState();
  return LuaCall<double>(L, lua_function_name_, mat_id, xyz);
}

} // namespace opensnlua
