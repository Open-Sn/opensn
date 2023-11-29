#include "lua/framework/math/functions/lua_scalar_spatial_function.h"
#include "framework/lua.h"
#include "framework/runtime.h"
#include "framework/console/console.h"
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
    lua_function_name_(params.GetParamValue<std::string>("lua_function_name"))
{
}

double
LuaScalarSpatialFunction::Evaluate(const opensn::Vector3& xyz) const
{
  const std::string fname = "LuaScalarSpatialFunction::Evaluate";

  lua_State* L = console.GetConsoleState();
  // Load lua function
  lua_getglobal(L, lua_function_name_.c_str());

  // Error check lua function
  if (not lua_isfunction(L, -1))
    ChiLogicalError("Attempted to access lua-function, " + lua_function_name_ +
                    ", but it seems the function could not be retrieved.");

  // Push arguments
  lua_pushnumber(L, xyz.x);
  lua_pushnumber(L, xyz.y);
  lua_pushnumber(L, xyz.z);

  // Call lua function
  // 3 arguments, 1 result (double), 0=original error object
  double lua_return;
  if (lua_pcall(L, 3, 1, 0) == 0)
  {
    LuaCheckNumberValue(fname, L, -1);
    lua_return = lua_tonumber(L, -1);
  }
  else
    ChiLogicalError("Attempted to call lua-function, " + lua_function_name_ +
                    ", but the call failed." + xyz.PrintStr());

  lua_pop(L, 1); // pop the double, or error code

  return lua_return;
}

} // namespace opensnlua
