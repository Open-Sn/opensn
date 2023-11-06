#ifdef OPENSN_WITH_LUA
#include "modules/LinearBoltzmannSolvers/C_DiscreteOrdinatesAdjointSolver/ResponseFunction/lbs_adj_response_function.h"

#include "framework/chi_lua.h"

#include "framework/chi_runtime.h"
#include "framework/console/chi_console.h"
#include "framework/logging/chi_log.h"

// ###################################################################
/** Calls the lua function associated with the response function and
 * returns a multigroup vector of the source values.*/
std::vector<double>
lbs::ResponseFunctionDesignation::GetMGResponse(const chi_mesh::Cell& cell,
                                                const size_t num_groups) const
{
  const std::string fname = __FUNCTION__;

  std::vector<double> response(num_groups, 0.0);

  //======================================== Utility lambdas
  auto PushVector3AsTable = [](lua_State* L, const chi_mesh::Vector3& vec)
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

  //============================================= Check response function given
  // Return default if none provided
  if (lua_functional.empty())
  {
    response.assign(num_groups, 1.0);
    return response;
  }

  //============================================= Load lua function
  lua_State* L = Chi::console.GetConsoleState();
  lua_getglobal(L, lua_functional.c_str());

  //============================================= Error check lua function
  if (not lua_isfunction(L, -1))
    throw std::logic_error(fname + " attempted to access lua-function, " + lua_functional +
                           ", but it seems the function"
                           " could not be retrieved.");

  //============================================= Push arguments
  PushVector3AsTable(L, cell.centroid_);
  lua_pushinteger(L, cell.material_id_); // 4 arguments on stack

  //============================================= Call lua function
  // 2 arguments, 1 result (table), 0=original error object
  std::vector<double> lua_return;
  if (lua_pcall(L, 2, 1, 0) == 0)
  {
    LuaCheckTableValue(fname, L, -1);
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
    throw std::logic_error(fname + " attempted to call lua-function, " + lua_functional +
                           ", but the call failed.");

  lua_pop(L, 1); // pop the table, or error code

  //============================================= Check return value
  if (lua_return.size() > response.size())
    throw std::logic_error(fname + " Call lua-function, " + lua_functional +
                           ", returned a vector of size " + std::to_string(response.size()) +
                           " which is greater than the number of groups " +
                           std::to_string(num_groups) + ".");

  for (size_t g = 0; g < num_groups; ++g)
    response[g] = lua_return[g];

  return response;
}
#endif
