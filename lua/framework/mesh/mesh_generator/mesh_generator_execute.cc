// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "lua/framework/lua.h"
#include "lua/framework/console/console.h"
#include "framework/mesh/mesh_generator/mesh_generator.h"
#include "framework/runtime.h"

using namespace opensn;

namespace opensnlua
{

int MeshGeneratorExecute(lua_State* L);

RegisterLuaFunctionInNamespace(MeshGeneratorExecute, mesh::MeshGenerator, Execute);

int
MeshGeneratorExecute(lua_State* L)
{
  const std::string fname = "mesh.MeshGenerator.Execute";
  LuaCheckArgs<size_t>(L, fname);

  auto handle = LuaArg<size_t>(L, 1);
  auto& generator = opensn::GetStackItem<MeshGenerator>(opensn::object_stack, handle, fname);
  generator.Execute();

  return LuaReturn(L);
}

} // namespace opensnlua
