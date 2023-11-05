#include "chi_lua.h"
#include <iostream>
#include "../surfacemesher.h"

#include "../../MeshHandler/chi_meshhandler.h"
#include "chi_runtime.h"
#include "chi_log.h"
#include "surfmesher_lua.h"
#include "console/chi_console.h"

RegisterLuaFunctionAsIs(chiSurfaceMesherSetProperty);
RegisterLuaConstantAsIs(MAX_AREA, chi_data_types::Varying(1));

//#############################################################################
/** Sets a property of a surface mesher.

\param PropertyNumber int Handle of the property to be set.
\param PropertyValue varying Value of the property.

Properties:\n
 MAX_AREA = Area constraint.\n

\ingroup LuaSurfaceMesher
\author Jan*/
int
chiSurfaceMesherSetProperty(lua_State* L)
{
  auto& cur_hndlr = chi_mesh::GetCurrentHandler();

  auto surf_mesher = cur_hndlr.GetSurfaceMesher();

  //================================================== Get property number
  int property_num = lua_tonumber(L, 1);

  //================================================== Area constraint
  if (property_num == 1) // MAX_AREA
  {
    Chi::log.Log0Warning() << "Deprecated and removed feature"
                              "property MAX_AREA in call"
                              " to chiSurfaceMesherSetProperty";
  }

  return 0;
}