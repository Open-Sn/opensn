#ifdef OPENSN_WITH_LUA
#include "framework/chi_lua.h"

#include "framework/ChiObjectFactory.h"
#include "framework/console/chi_console.h"

namespace chi::lua_utils
{

int chiMakeObject(lua_State* L);
int chiMakeObjectType(lua_State* L);

RegisterLuaFunctionAsIs(chiMakeObject);
RegisterLuaFunctionAsIs(chiMakeObjectType);

// #############################################################################
/**Generic lua routine for the creation of objects.
 * \param params ParameterBlock A single block tree that requires a parameter
 *  called chi_obj_type that indicates the type of object to make.

*/
int
chiMakeObject(lua_State* L)
{
  const std::string fname = __FUNCTION__;
  const int num_args = lua_gettop(L);
  if (num_args != 1) LuaPostArgAmountError(fname, 2, num_args);

  LuaCheckTableValue(fname, L, 1);

  const auto params = chi_lua::TableParserAsParameterBlock::ParseTable(L, 1);

  const auto& object_maker = ChiObjectFactory::GetInstance();
  const size_t handle = object_maker.MakeRegisteredObject(params);

  const std::string type = params.GetParamValue<std::string>("chi_obj_type");

  lua_pushinteger(L, static_cast<lua_Integer>(handle));
  return 1;
}

// #############################################################################
/**Generic lua routine for the creation of objects.
 * \param type string The type to create.
 * \param params ParameterBlock A single block tree.

*/
int
chiMakeObjectType(lua_State* L)
{
  const std::string fname = __FUNCTION__;
  const int num_args = lua_gettop(L);
  if (num_args != 2) LuaPostArgAmountError(fname, 2, num_args);

  LuaCheckStringValue(fname, L, 1);
  LuaCheckTableValue(fname, L, 2);

  const std::string type = lua_tostring(L, 1);
  const auto params = chi_lua::TableParserAsParameterBlock::ParseTable(L, 2);

  const auto& object_maker = ChiObjectFactory::GetInstance();
  const size_t handle = object_maker.MakeRegisteredObjectOfType(type, params);

  lua_pushinteger(L, static_cast<lua_Integer>(handle));
  return 1;
}

} // namespace chi::lua_utils
#endif
