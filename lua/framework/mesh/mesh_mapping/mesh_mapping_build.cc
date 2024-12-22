// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "lua/framework/lua.h"
#include "lua/framework/console/console.h"
#include "framework/mesh/mesh_generator/mesh_generator.h"
#include "framework/runtime.h"
#include "framework/mesh/mesh_mapping/mesh_mapping.h"

using namespace opensn;

namespace opensnlua
{

int MeshMappingBuild(lua_State* L);

RegisterLuaFunctionInNamespace(MeshMappingBuild, mesh::MeshMapping, Build);

int
MeshMappingBuild(lua_State* L)
{
  // This mesh search is really bad, but it's good enough for now. Once
  // #435 gets in, we can take the mesh generator handles directly.
  // It assumes that the fine mesh is created first and the coarse
  // mesh is generated second.
  if (mesh_stack.size() != 2)
    throw std::logic_error(
      "In mesh.MeshMapping.Build: Two meshes are not available. This build assumes that the first "
      "mesh in the stack is the fine mesh and the second mesh in the stack is the coarse mesh.");

  MeshMapping mesh_mapping;
  mesh_mapping.Build(*mesh_stack[0], *mesh_stack[1]);

  return LuaReturn(L);
}

} // namespace opensnlua
