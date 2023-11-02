#include "chi_lua.h"

#include "mesh/VolumeMesher/chi_volumemesher.h"
#include "volumemesher_lua.h"
#include "console/chi_console.h"

RegisterLuaFunctionAsIs(chiVolumeMesherSetupOrthogonalBoundaries);

RegisterLuaConstant(OrthoBoundaryID, XMAX,chi_data_types::Varying(0));
RegisterLuaConstant(OrthoBoundaryID, XMIN,chi_data_types::Varying(1));
RegisterLuaConstant(OrthoBoundaryID, YMAX,chi_data_types::Varying(2));
RegisterLuaConstant(OrthoBoundaryID, YMIN,chi_data_types::Varying(3));
RegisterLuaConstant(OrthoBoundaryID, ZMAX,chi_data_types::Varying(4));
RegisterLuaConstant(OrthoBoundaryID, ZMIN,chi_data_types::Varying(5));

//###################################################################
/** Sets boundary numbers on boundaries orthogonal to the cardinal directions
 * as xmax=0, xmin=1, ymax=2, ymin=3, zmax=4, zmin=5.
 *
\ingroup LuaVolumeMesher
 */
int chiVolumeMesherSetupOrthogonalBoundaries(lua_State* L)
{
  chi_mesh::VolumeMesher::SetupOrthogonalBoundaries();
  return 0;
}