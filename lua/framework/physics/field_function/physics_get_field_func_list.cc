#include "framework/lua.h"
#include "framework/runtime.h"
#include "framework/physics/solver_base/solver.h"
#include "framework/field_functions/field_function_grid_based.h"
#include "framework/logging/log.h"
#include "field_functions_lua.h"
#include "framework/console/console.h"

using namespace opensn;

RegisterLuaFunctionAsIs(GetFieldFunctionHandleByName);

int
GetFieldFunctionHandleByName(lua_State* L)
{
  const std::string fname = __FUNCTION__;
  int num_args = lua_gettop(L);
  if (num_args != 1)
    LuaPostArgAmountError(fname, 1, num_args);

  LuaCheckNilValue(fname, L, 1);
  LuaCheckStringValue(fname, L, 1);

  const std::string ff_name = lua_tostring(L, 1);

  size_t ff_handle_counter = 0;
  std::vector<size_t> handles_that_matched;
  for (const auto& pff : opensn::field_function_stack)
  {
    if (pff->TextName() == ff_name)
      handles_that_matched.emplace_back(ff_handle_counter);
    ++ff_handle_counter;
  }

  size_t num_handles = handles_that_matched.size();

  if (num_handles == 0)
  {
    opensn::log.Log0Warning() << fname << ": No field-functions were found that "
                              << "matched the requested name:\"" << ff_name
                              << "\". A null handle will "
                              << "be returned." << std::endl;

    return 0;
  }

  if (num_handles > 1)
    opensn::log.Log0Warning() << fname << ": A total of " << num_handles
                              << " field-functions were found that matched the "
                              << " requested name. Only the first match will be "
                              << " returned.";

  lua_pushinteger(L, static_cast<lua_Integer>(handles_that_matched.front()));
  return 1;
}
