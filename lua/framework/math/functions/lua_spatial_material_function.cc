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
  // Utility lambdas
  auto PushVector3AsTable = [](lua_State* L, const Vector3& vec)
  {
    lua_newtable(L);

    lua_pushstring(L, "x");
    lua_pushnumber(L, vec.x);
    lua_settable(L, -3);

    lua_pushstring(L, "y");
    lua_pushnumber(L, vec.y);
    lua_settable(L, -3);

    lua_pushstring(L, "z");
    lua_pushnumber(L, vec.z);
    lua_settable(L, -3);
  };

  // Check response function given
  // Return default if none provided
  if (lua_function_name_.empty())
    return std::vector<double>(num_components, 1.0);

  // Load lua function
  lua_State* L = console.GetConsoleState();
  lua_getglobal(L, lua_function_name_.c_str());

  // Error check lua function
  if (not lua_isfunction(L, -1))
    OpenSnLogicalError("Attempted to access lua-function, " + lua_function_name_ +
                       ", but it seems the function could not be retrieved.");

  // Push arguments
  PushVector3AsTable(L, xyz);
  lua_pushinteger(L, mat_id); // 4 arguments on stack

  // Call lua function
  // 2 arguments, 1 result (table), 0=original error object
  std::vector<double> lua_return;
  if (lua_pcall(L, 2, 1, 0) == 0)
  {
    LuaCheckTableValue(__FUNCTION__, L, -1);
    const size_t table_length = lua_rawlen(L, -1);
    lua_return.reserve(table_length);
    for (size_t i = 0; i < table_length; ++i)
    {
      lua_pushinteger(L, static_cast<lua_Integer>(i) + 1);
      lua_gettable(L, -2);
      lua_return.push_back(lua_tonumber(L, -1));
      lua_pop(L, 1);
    }
  }
  else
    OpenSnLogicalError("attempted to call lua-function, " + lua_function_name_ +
                       ", but the call failed.");

  lua_pop(L, 1); // pop the table, or error code

  // Check return value
  OpenSnLogicalErrorIf(lua_return.size() != num_components,
                       "Call to lua function " + lua_function_name_ + " returned a vector of " +
                         "size " + std::to_string(lua_return.size()) +
                         ", which is not the same as " + "the number of groups " +
                         std::to_string(num_components) + ".");

  return lua_return;
}

} // namespace opensnlua
