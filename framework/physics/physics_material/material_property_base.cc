#include "framework/physics/physics_material/material_property_base.h"

#ifdef OPENSN_WITH_LUA

namespace opensn
{

void
MaterialProperty::PushLuaTable(lua_State* L) const
{
  lua_newtable(L);
  lua_pushstring(L, "is_empty");
  lua_pushboolean(L, true);
  lua_settable(L, -3);
}

} // namespace opensn
#endif
