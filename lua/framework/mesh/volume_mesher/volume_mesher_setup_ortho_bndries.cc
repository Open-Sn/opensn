#include "framework/lua.h"

#include "framework/mesh/volume_mesher/volume_mesher.h"
#include "volume_mesher_lua.h"
#include "framework/console/console.h"

RegisterLuaFunctionAsIs(chiVolumeMesherSetupOrthogonalBoundaries);

RegisterLuaConstant(OrthoBoundaryID, XMAX, chi_data_types::Varying(0));
RegisterLuaConstant(OrthoBoundaryID, XMIN, chi_data_types::Varying(1));
RegisterLuaConstant(OrthoBoundaryID, YMAX, chi_data_types::Varying(2));
RegisterLuaConstant(OrthoBoundaryID, YMIN, chi_data_types::Varying(3));
RegisterLuaConstant(OrthoBoundaryID, ZMAX, chi_data_types::Varying(4));
RegisterLuaConstant(OrthoBoundaryID, ZMIN, chi_data_types::Varying(5));

int
chiVolumeMesherSetupOrthogonalBoundaries(lua_State* L)
{
  chi_mesh::VolumeMesher::SetupOrthogonalBoundaries();
  return 0;
}
