// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "lbs_common_lua_functions.h"
#include "lua/framework/lua.h"

#define RegisterTable(x)                                                                           \
  lua_newtable(L);                                                                                 \
  lua_setglobal(L, #x)

#define RegisterNumberValueToTable(const_name, const_value, namespace_name)                        \
  lua_getglobal(L, #namespace_name);                                                               \
  lua_pushstring(L, #const_name);                                                                  \
  lua_pushnumber(L, const_value);                                                                  \
  lua_settable(L, -3);                                                                             \
  lua_pop(L, 1)

using namespace opensn;

namespace opensnlua::lbs
{
void
RegisterLuaEntities(lua_State* L)
{
  LuaSetGlobal(L, "DISCRETIZATION_METHOD", 1);
  LuaSetGlobal(L, "PWLD", 3);
  LuaSetGlobal(L, "PWLD1D", 4);
  LuaSetGlobal(L, "PWLD2D", 5);
  LuaSetGlobal(L, "PWLD3D", 6);
  LuaSetGlobal(L, "PARTITION_METHOD", 2);
  LuaSetGlobal(L, "SERIAL", 1);
  LuaSetGlobal(L, "FROM_SURFACE", 2);
  LuaSetGlobal(L, "BOUNDARY_CONDITION", 3);
  LuaSetGlobal(L, "XMAX", 31);
  LuaSetGlobal(L, "XMIN", 32);
  LuaSetGlobal(L, "YMAX", 33);
  LuaSetGlobal(L, "YMIN", 34);
  LuaSetGlobal(L, "ZMAX", 35);
  LuaSetGlobal(L, "ZMIN", 36);
  LuaSetGlobal(L, "SCATTERING_ORDER", 4);
  LuaSetGlobal(L, "MAX_MPI_MESSAGE_SIZE", 5);
  LuaSetGlobal(L, "READ_RESTART_DATA", 6);
  LuaSetGlobal(L, "WRITE_RESTART_DATA", 7);
  LuaSetGlobal(L, "SAVE_ANGULAR_FLUX", 8);
  LuaSetGlobal(L, "USE_SOURCE_MOMENTS", 9);
  LuaSetGlobal(L, "VERBOSE_INNER_ITERATIONS", 10);
  LuaSetGlobal(L, "VERBOSE_OUTER_ITERATIONS", 11);
  LuaSetGlobal(L, "USE_PRECURSORS", 12);

  RegisterTable(LBSBoundaryTypes);
  RegisterNumberValueToTable(VACUUM, 1, LBSBoundaryTypes);
  RegisterNumberValueToTable(INCIDENT_ISOTROPIC, 2, LBSBoundaryTypes);
  RegisterNumberValueToTable(REFLECTING, 3, LBSBoundaryTypes);
  RegisterNumberValueToTable(INCIDENT_ANISTROPIC_HETEROGENEOUS, 4, LBSBoundaryTypes);

  RegisterTable(LBSGroupset);
  RegisterNumberValueToTable(ANGLE_AGG_SINGLE, 1, LBSGroupset);
  RegisterNumberValueToTable(ANGLE_AGG_POLAR, 2, LBSGroupset);
  RegisterNumberValueToTable(ANGLE_AGG_AZIMUTHAL, 3, LBSGroupset);
  LuaSetGlobal(L, "KRYLOV_RICHARDSON", 5);
  LuaSetGlobal(L, "KRYLOV_RICHARDSON_CYCLES", 6);
  LuaSetGlobal(L, "KRYLOV_GMRES", 7);
  LuaSetGlobal(L, "KRYLOV_GMRES_CYCLES", 8);
  LuaSetGlobal(L, "KRYLOV_BICGSTAB", 9);
  LuaSetGlobal(L, "KRYLOV_BICGSTAB_CYCLES", 10);
}

} // namespace opensnlua::lbs
