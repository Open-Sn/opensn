#include "framework/lua.h"
#include <iostream>
#include "framework/physics/physics_material/physics_material.h"

#include "framework/runtime.h"
#include "framework/logging/log.h"
#include "physics_lua_utils.h"
#include "framework/console/console.h"

using namespace opensn;

RegisterLuaFunctionAsIs(chiPhysicsAddMaterial);

int
chiPhysicsAddMaterial(lua_State* L)
{
  int numArgs = lua_gettop(L);

  auto new_material = std::make_shared<Material>();
  if (numArgs == 1)
  {
    const char* temp = lua_tostring(L, 1);
    new_material->name_ = std::string(temp);
  }

  opensn::Chi::material_stack.push_back(new_material);

  const size_t index = opensn::Chi::material_stack.size() - 1;
  lua_pushnumber(L, static_cast<lua_Number>(index));

  opensn::Chi::log.Log0Verbose1() << "New material added at index " << index << " with name \""
                                  << new_material->name_ << "\"";

  return 1;
}
