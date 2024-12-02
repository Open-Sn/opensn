// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "lua/lib/math/functions/lua_scalar_material_function.h"
#include "framework/runtime.h"
#include "lua/lib/console.h"
#include "lua/lib/types.h"
#include "LuaBridge/LuaBridge.h"
#include <lua.h>

using namespace opensn;

namespace opensnlua
{

OpenSnRegisterObjectInNamespace(opensn, LuaScalarMaterialFunction);

InputParameters
LuaScalarMaterialFunction::GetInputParameters()
{
  InputParameters params = ScalarMaterialFunction::GetInputParameters();
  params.AddRequiredParameter<std::string>("function_name", "Name of the lua function");
  return params;
}

std::shared_ptr<LuaScalarMaterialFunction>
LuaScalarMaterialFunction::Create(const opensn::ParameterBlock& params)
{
  auto& factory = opensn::ObjectFactory::GetInstance();
  return factory.Create<LuaScalarMaterialFunction>("opensn::LuaScalarMaterialFunction", params);
}

LuaScalarMaterialFunction::LuaScalarMaterialFunction(const InputParameters& params)
  : ScalarMaterialFunction(params),
    function_name_(params.GetParamValue<std::string>("function_name"))
{
}

double
LuaScalarMaterialFunction::Evaluate(double val, int mat_id) const
{
  lua_State* L = console.GetConsoleState();
  auto fn = luabridge::getGlobal(L, function_name_.c_str());
  return fn(val, mat_id)[0];
}

} // namespace opensnlua
