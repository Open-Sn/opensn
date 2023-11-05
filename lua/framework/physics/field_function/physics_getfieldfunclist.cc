#include "opensn/framework/chi_lua.h"

#include "opensn/framework/chi_runtime.h"

#include "opensn/framework/physics/SolverBase/chi_solver.h"
#include "opensn/framework/physics/FieldFunction/fieldfunction_gridbased.h"

#include "opensn/framework/logging/chi_log.h"
#include "fieldfunctions_lua.h"
#include "opensn/framework/console/chi_console.h"

RegisterLuaFunctionAsIs(chiGetFieldFunctionHandleByName);

//###################################################################
/**Gets a field-function handle by name.
\param FFname string Name of the field function.

\return handle If the field-function was found and a handle identified the valid
               handle will be returned (i.e., a natural number >= 0). If the
               field-function by the given name was not found then the function
               will return null.

\ingroup LuaFieldFunc
 */
int
chiGetFieldFunctionHandleByName(lua_State* L)
{
  const std::string fname = __FUNCTION__;
  int num_args = lua_gettop(L);
  if (num_args != 1) LuaPostArgAmountError(fname, 1, num_args);

  LuaCheckNilValue(fname, L, 1);
  LuaCheckStringValue(fname, L, 1);

  const std::string ff_name = lua_tostring(L, 1);

  size_t ff_handle_counter = 0;
  std::vector<size_t> handles_that_matched;
  for (const auto& pff : Chi::field_function_stack)
  {
    if (pff->TextName() == ff_name) handles_that_matched.emplace_back(ff_handle_counter);
    ++ff_handle_counter;
  }

  size_t num_handles = handles_that_matched.size();

  if (num_handles == 0)
  {
    Chi::log.Log0Warning() << fname << ": No field-functions were found that "
                           << "matched the requested name:\"" << ff_name
                           << "\". A null handle will "
                           << "be returned." << std::endl;

    return 0;
  }

  if (num_handles > 1)
    Chi::log.Log0Warning() << fname << ": A total of " << num_handles
                           << " field-functions were found that matched the "
                           << " requested name. Only the first match will be "
                           << " returned.";

  lua_pushinteger(L, static_cast<lua_Integer>(handles_that_matched.front()));
  return 1;
}
