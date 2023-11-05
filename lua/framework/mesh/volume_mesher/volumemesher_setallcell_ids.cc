#include "opensn/framework/chi_lua.h"

#include "opensn/framework/mesh/VolumeMesher/chi_volumemesher.h"
#include "volumemesher_lua.h"
#include "opensn/framework/console/chi_console.h"

RegisterLuaFunctionAsIs(chiVolumeMesherSetMatIDToAll);

//###################################################################
/** Sets all cell-material id's to the supplied value.
\param material_id int The id.
\ingroup LuaVolumeMesher
*/
int
chiVolumeMesherSetMatIDToAll(lua_State* L)
{
  int num_args = lua_gettop(L);
  if (num_args != 1) LuaPostArgAmountError(__FUNCTION__, 1, num_args);

  LuaCheckNilValue(__FUNCTION__, L, 1);

  int mat_id = lua_tonumber(L, 1);

  chi_mesh::VolumeMesher::SetMatIDToAll(mat_id);
  return 0;
}
