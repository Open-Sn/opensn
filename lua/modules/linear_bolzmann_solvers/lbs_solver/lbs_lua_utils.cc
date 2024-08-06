// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "lua/framework/lua.h"

using namespace opensn;

namespace opensnlua
{

void
RegisterLuaEntitiesLBS(lua_State* L)
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
  LuaSetGlobal(L, "XMIN", 31);
  LuaSetGlobal(L, "XMAX", 32);
  LuaSetGlobal(L, "YMIN", 33);
  LuaSetGlobal(L, "YMAX", 34);
  LuaSetGlobal(L, "ZMIN", 35);
  LuaSetGlobal(L, "ZMAX", 36);
  LuaSetGlobal(L, "SCATTERING_ORDER", 4);
  LuaSetGlobal(L, "MAX_MPI_MESSAGE_SIZE", 5);
  LuaSetGlobal(L, "READ_RESTART_DATA", 6);
  LuaSetGlobal(L, "WRITE_RESTART_DATA", 7);
  LuaSetGlobal(L, "SAVE_ANGULAR_FLUX", 8);
  LuaSetGlobal(L, "USE_SOURCE_MOMENTS", 9);
  LuaSetGlobal(L, "VERBOSE_INNER_ITERATIONS", 10);
  LuaSetGlobal(L, "VERBOSE_OUTER_ITERATIONS", 11);
  LuaSetGlobal(L, "USE_PRECURSORS", 12);

  LuaRegisterTable(L, "LBSBoundaryTypes");
  LuaRegisterTableConstant(L, "LBSBoundaryTypes", "VACUUM", 1);
  LuaRegisterTableConstant(L, "LBSBoundaryTypes", "INCIDENT_ISOTROPIC", 2);
  LuaRegisterTableConstant(L, "LBSBoundaryTypes", "REFLECTING", 3);
  LuaRegisterTableConstant(L, "LBSBoundaryTypes", "INCIDENT_ANISTROPIC_HETEROGENEOUS", 4);

  LuaRegisterTable(L, "LBSGroupset");
  LuaRegisterTableConstant(L, "LBSGroupset", "ANGLE_AGG_SINGLE", 1);
  LuaRegisterTableConstant(L, "LBSGroupset", "ANGLE_AGG_POLAR", 2);
  LuaRegisterTableConstant(L, "LBSGroupset", "ANGLE_AGG_AZIMUTHAL", 3);
  LuaSetGlobal(L, "KRYLOV_RICHARDSON", 5);
  LuaSetGlobal(L, "KRYLOV_RICHARDSON_CYCLES", 6);
  LuaSetGlobal(L, "KRYLOV_GMRES", 7);
  LuaSetGlobal(L, "KRYLOV_GMRES_CYCLES", 8);
  LuaSetGlobal(L, "KRYLOV_BICGSTAB", 9);
  LuaSetGlobal(L, "KRYLOV_BICGSTAB_CYCLES", 10);
}

} // namespace opensnlua
