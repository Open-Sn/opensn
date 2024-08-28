// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "lua/framework/math/functions/lua_scalar_material_function.h"
#include "lua/framework/lua.h"
#include "lua/framework/console/console.h"
#include "framework/runtime.h"
#include "framework/object_factory.h"

using namespace opensn;

namespace opensnlua
{

OpenSnRegisterObjectInNamespace(opensn, LuaScalarMaterialFunction);

InputParameters
LuaScalarMaterialFunction::GetInputParameters()
{
  InputParameters params = ScalarMaterialFunction::GetInputParameters();
  params.AddRequiredParameter<std::string>("lua_function_name", "Name of the lua function");
  return params;
}

LuaScalarMaterialFunction::LuaScalarMaterialFunction(const InputParameters& params)
  : ScalarMaterialFunction(params),
    lua_function_name_(params.GetParamValue<std::string>("lua_function_name"))
{
}

double
LuaScalarMaterialFunction::Evaluate(double val, int mat_id) const
{
  lua_State* L = console.GetConsoleState();
  return LuaCall<double>(L, lua_function_name_, val, mat_id);
}

} // namespace opensnlua
