#pragma once

#include "opensn/framework/physics/PhysicsMaterial/material_property_base.h"

namespace chi_physics
{

//###################################################################
/**Simple scalar material property.*/
class ScalarValue : public chi_physics::MaterialProperty
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

} // namespace chi_physics
