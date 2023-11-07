#include "framework/lua.h"

#include "framework/mesh/mesh_cutting/mesh_cutting.h"
#include "framework/mesh/mesh.h"
#include "framework/mesh/mesh_handler/mesh_handler.h"

#include "mesh_cutting_lua.h"
#include "framework/console/console.h"

RegisterLuaFunctionAsIs(chiCutMesh);

int
chiCutMesh(lua_State* L)
{
  auto fname = __FUNCTION__;
  int num_args = lua_gettop(L);
  if (num_args < 2) LuaPostArgAmountError(fname, 2, num_args);

  LuaCheckNilValue(fname, L, 1);
  LuaCheckNilValue(fname, L, 2);
  if (num_args >= 3) LuaCheckNilValue(fname, L, 3);
  if (num_args == 4) LuaCheckNilValue(fname, L, 4);

  std::vector<double> p_raw;
  std::vector<double> n_raw;
  double merge_tolerance = 1.0e-3;
  double float_compare = 1.0e-10;

  LuaPopulateVectorFrom1DArray(fname, L, 1, p_raw);
  LuaPopulateVectorFrom1DArray(fname, L, 2, n_raw);
  if (num_args >= 3) merge_tolerance = lua_tonumber(L, 3);
  if (num_args == 4) float_compare = lua_tonumber(L, 4);

  auto& handler = chi_mesh::GetCurrentHandler();

  auto& grid = handler.GetGrid();

  chi_mesh::Vector3 p(p_raw[0], p_raw[1], p_raw[2]);
  chi_mesh::Vector3 n(n_raw[0], n_raw[1], n_raw[2]);

  chi_mesh::mesh_cutting::CutMeshWithPlane(*grid, p, n, merge_tolerance, float_compare);

  return 0;
}
