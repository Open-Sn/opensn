// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "lua/framework/mesh/mesh_compute_vol.h"
#include "framework/mesh/mesh_continuum/mesh_continuum.h"
#include "framework/runtime.h"
#include "framework/logging/log.h"
#include "lua/framework/console/console.h"

namespace opensnlua
{

RegisterLuaFunctionInNamespace(MeshComputeVolumePerMaterialID, mesh, ComputeVolumePerMaterialID);

using namespace opensn;

int
MeshComputeVolumePerMaterialID(lua_State* L)
{
  auto curr_mesh = GetCurrentMesh();
  curr_mesh->ComputeVolumePerMaterialID();

  return LuaReturn(L);
}

} // namespace opensnlua
