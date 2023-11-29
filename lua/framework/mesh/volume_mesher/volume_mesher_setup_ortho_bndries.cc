#include "framework/lua.h"

#include "framework/mesh/volume_mesher/volume_mesher.h"
#include "volume_mesher_lua.h"
#include "framework/console/console.h"

using namespace opensn;

RegisterLuaFunctionAsIs(VolumeMesherSetupOrthogonalBoundaries);

RegisterLuaConstant(OrthoBoundaryID, XMAX, Varying(0));
RegisterLuaConstant(OrthoBoundaryID, XMIN, Varying(1));
RegisterLuaConstant(OrthoBoundaryID, YMAX, Varying(2));
RegisterLuaConstant(OrthoBoundaryID, YMIN, Varying(3));
RegisterLuaConstant(OrthoBoundaryID, ZMAX, Varying(4));
RegisterLuaConstant(OrthoBoundaryID, ZMIN, Varying(5));

int
VolumeMesherSetupOrthogonalBoundaries(lua_State* L)
{
  VolumeMesher::SetupOrthogonalBoundaries();
  return 0;
}
