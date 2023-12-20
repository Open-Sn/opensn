#pragma once

#include "framework/physics/physics_material/material_property_base.h"

namespace opensn
{

/**Simple scalar material property.*/
class ScalarValue : public MaterialProperty
{
public:
  double value_ = 1.0;

  ScalarValue() : MaterialProperty(PropertyType::SCALAR_VALUE) {}

  double GetScalarValue() override { return value_; }
#ifdef OPENSN_WITH_LUA
  void PushLuaTable(lua_State* L) const override
  {
    lua_newtable(L);
    lua_pushstring(L, "is_empty");
    lua_pushboolean(L, false);
    lua_settable(L, -3);

    lua_pushstring(L, "value");
    lua_pushnumber(L, value_);
    lua_settable(L, -3);
  }
#endif
};

} // namespace opensn
