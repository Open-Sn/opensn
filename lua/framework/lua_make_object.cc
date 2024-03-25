#include "lua/framework/lua.h"

#include "framework/object_factory.h"
#include "lua/framework/console/console.h"

namespace opensnlua
{

/**Generic lua routine for the creation of objects.
 * \param params ParameterBlock A single block tree that requires a parameter
 *  called obj_type that indicates the type of object to make.
 */
int MakeObject(lua_State* L);

/**Generic lua routine for the creation of objects.
 * \param type string The type to create.
 * \param params ParameterBlock A single block tree.
 */
int MakeObjectType(lua_State* L);

RegisterLuaFunctionAsIs(MakeObject);
RegisterLuaFunctionAsIs(MakeObjectType);

int
MakeObject(lua_State* L)
{
  const std::string fname = __FUNCTION__;
  const int num_args = lua_gettop(L);
  if (num_args != 1)
    LuaPostArgAmountError(fname, 2, num_args);

  LuaCheckTableValue(fname, L, 1);

  const auto params = TableParserAsParameterBlock::ParseTable(L, 1);

  const auto& object_maker = opensn::ObjectFactory::GetInstance();
  const auto handle = object_maker.MakeRegisteredObject(params);
  return LuaReturn(L, handle);
}

int
MakeObjectType(lua_State* L)
{
  const std::string fname = __FUNCTION__;
  const int num_args = lua_gettop(L);
  if (num_args != 2)
    LuaPostArgAmountError(fname, 2, num_args);

  LuaCheckTableValue(fname, L, 2);

  const auto type = LuaArg<std::string>(L, 1);
  const auto params = TableParserAsParameterBlock::ParseTable(L, 2);

  const auto& object_maker = opensn::ObjectFactory::GetInstance();
  const auto handle = object_maker.MakeRegisteredObjectOfType(type, params);
  return LuaReturn(L, handle);
}

} // namespace opensnlua
