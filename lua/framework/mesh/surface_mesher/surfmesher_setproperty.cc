#include "framework/lua.h"
#include <iostream>
#include "framework/mesh/surface_mesher/surface_mesher.h"

#include "framework/mesh/mesh_handler/mesh_handler.h"
#include "framework/runtime.h"
#include "framework/logging/log.h"
#include "surfmesher_lua.h"
#include "framework/console/console.h"

RegisterLuaFunctionAsIs(chiSurfaceMesherSetProperty);
RegisterLuaConstantAsIs(MAX_AREA, chi_data_types::Varying(1));

int
chiSurfaceMesherSetProperty(lua_State* L)
{
  auto& cur_hndlr = chi_mesh::GetCurrentHandler();

  auto surf_mesher = cur_hndlr.GetSurfaceMesher();

  // Get property number
  int property_num = lua_tonumber(L, 1);

  // Area constraint
  if (property_num == 1) // MAX_AREA
  {
    Chi::log.Log0Warning() << "Deprecated and removed feature"
                              "property MAX_AREA in call"
                              " to chiSurfaceMesherSetProperty";
  }

  return 0;
}
