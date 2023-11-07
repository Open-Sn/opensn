#pragma once

#include "framework/physics/PhysicsMaterial/material_property_base.h"

namespace chi_physics
{

/** Basic thermal conductivity material property.*/
class IsotropicMultiGrpSource : public chi_physics::MaterialProperty
{
public:
  std::vector<double> source_value_g_;

  IsotropicMultiGrpSource() : MaterialProperty(PropertyType::ISOTROPIC_MG_SOURCE) {}

#ifdef OPENSN_WITH_LUA
  void PushLuaTable(lua_State* L) const override
  {
    lua_newtable(L);
    lua_pushstring(L, "is_empty");
    lua_pushboolean(L, false);
    lua_settable(L, -3);

    lua_pushstring(L, "G");
    lua_pushnumber(L, source_value_g_.size());
    lua_settable(L, -3);

    lua_pushstring(L, "source_value_g");
    lua_newtable(L);
    int g = 0;
    for (auto val : source_value_g_)
    {
      ++g;
      lua_pushnumber(L, g);
      lua_pushnumber(L, val);
      lua_settable(L, -3);
    }
    lua_settable(L, -3);
  }
#endif
};

} // namespace chi_physics
