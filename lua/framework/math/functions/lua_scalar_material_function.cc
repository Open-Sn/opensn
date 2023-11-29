#include "lua/framework/math/functions/lua_scalar_material_function.h"
#include "framework/lua.h"
#include "framework/runtime.h"
#include "framework/console/console.h"
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
  double ret_val = 0.0;

  lua_getglobal(L, lua_function_name_.c_str());
  lua_pushnumber(L, val);
  lua_pushnumber(L, mat_id);

  // 2 arguments, 1 result, 0=original error object
  if (lua_pcall(L, 2, 1, 0) == 0) { ret_val = lua_tonumber(L, -1); }
  lua_pop(L, 1);

  return ret_val;
}

} // namespace opensnlua
