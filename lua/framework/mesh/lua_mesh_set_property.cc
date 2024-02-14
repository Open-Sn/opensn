#include "framework/lua.h"
#include "framework/mesh/mesh_continuum/mesh_continuum.h"
#include "framework/utils/timer.h"
#include "framework/runtime.h"
#include "framework/mesh/logical_volume/logical_volume.h"
#include "lua/framework/mesh/lua_mesh_set_property.h"
#include "lua/framework/mesh/lua_mesh_set_ids.h"
#include "framework/console/console.h"
#include <iostream>

using namespace opensn;

RegisterLuaFunctionNamespace(MeshSetProperty, mesh, SetProperty);

enum MeshProperty
{
  MATID_FROM_LUA_FUNCTION = 13,
  BNDRYID_FROM_LUA_FUNCTION = 14
};

RegisterLuaConstantAsIs(MATID_FROM_LUA_FUNCTION, Varying(13));
RegisterLuaConstantAsIs(BNDRYID_FROM_LUA_FUNCTION, Varying(14));

int
MeshSetProperty(lua_State* L)
{
  const std::string fname = "mesh.SetProperty";

  // Get property index
  const int num_args = lua_gettop(L);
  if (num_args < 2)
    LuaPostArgAmountError(fname, 2, num_args);

  LuaCheckNilValue(fname, L, 1);
  LuaCheckNilValue(fname, L, 2);

  int property_index = lua_tonumber(L, 1);

  if (property_index == MeshProperty::MATID_FROM_LUA_FUNCTION)
  {
    LuaCheckStringValue(fname, L, 2);

    const std::string lua_fname = lua_tostring(L, 2);

    SetMatIDFromLuaFunction(lua_fname);
  }
  else if (property_index == MeshProperty::BNDRYID_FROM_LUA_FUNCTION)
  {
    LuaCheckStringValue(fname, L, 2);

    const std::string lua_fname = lua_tostring(L, 2);

    SetBoundaryIDFromLuaFunction(lua_fname);
  }
  else
  {
    opensn::log.LogAllError() << "Invalid property specified " << property_index
                              << " in call to VolumeMesherSetProperty().";
    opensn::Exit(EXIT_FAILURE);
  }

  return 0;
}
