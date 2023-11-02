#ifndef CHITECH_UNPARTITION_MESH_LUA_UTILS_H
#define CHITECH_UNPARTITION_MESH_LUA_UTILS_H

#include "mesh/UnpartitionedMesh/chi_unpartitioned_mesh.h"

#include "chi_lua.h"

namespace chi_mesh::unpartition_mesh_lua_utils
{
  //create.cc
  int chiCreateEmptyUnpartitionedMesh(lua_State* L);
  int chiDestroyUnpartitionedMesh(lua_State* L);

  int chiUnpartitionedMeshFromVTU(lua_State* L);
  int chiUnpartitionedMeshFromPVTU(lua_State* L);
  int chiUnpartitionedMeshFromEnsightGold(lua_State* L);
  int chiUnpartitionedMeshFromWavefrontOBJ(lua_State* L);
  int chiUnpartitionedMeshFromMshFormat(lua_State* L);
  int chiUnpartitionedMeshFromExodusII(lua_State* L);

  //basic_operations.cc
  int chiUnpartitionedMeshUploadVertex(lua_State* L);
  int chiUnpartitionedMeshUploadCell(lua_State* L);
  int chiUnpartitionedMeshFinalizeEmpty(lua_State* L);

}//namespace chi_mesh::unpartition_mesh_lua_utils



#endif //CHITECH_UNPARTITION_MESH_LUA_UTILS_H
