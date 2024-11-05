// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "lua/framework/lua.h"
#include "lua/framework/console/console.h"
#include "framework/object_factory.h"

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

RegisterLuaFunction(MakeObject);
RegisterLuaFunction(MakeObjectType);

int
MakeObject(lua_State* L)
{
  const std::string fname = __FUNCTION__;
  LuaCheckArgs<opensn::ParameterBlock>(L, fname);

  const auto params = LuaArg<opensn::ParameterBlock>(L, 1);
  const auto& object_maker = opensn::ObjectFactory::Instance();
  const auto handle = object_maker.MakeRegisteredObject(params);
  return LuaReturn(L, handle);
}

int
MakeObjectType(lua_State* L)
{
  const std::string fname = __FUNCTION__;
  LuaCheckArgs<std::string, opensn::ParameterBlock>(L, fname);

  const auto type = LuaArg<std::string>(L, 1);
  const auto params = LuaArg<opensn::ParameterBlock>(L, 2);
  const auto& object_maker = opensn::ObjectFactory::Instance();
  const auto handle = object_maker.MakeRegisteredObjectOfType(type, params);
  return LuaReturn(L, handle);
}

} // namespace opensnlua
