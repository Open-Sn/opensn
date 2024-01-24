#include "lbsadj_lua_utils.h"

#include "modules/linear_boltzmann_solvers/c_discrete_ordinates_adjoint_solver/lbs_adj_solver.h"
#include "lua/framework/math/functions/lua_spatial_material_function.h"
#include "framework/runtime.h"
#include "framework/mesh/logical_volume/logical_volume.h"

#include "framework/console/console.h"

using namespace opensn;

namespace opensnlua::lbs
{

namespace
{

std::shared_ptr<LuaSpatialMaterialFunction>
CreateResponseFunction(const std::string& function_name)
{
  ParameterBlock blk;
  blk.AddParameter("lua_function_name", function_name);
  InputParameters params = LuaSpatialMaterialFunction::GetInputParameters();
  params.AssignParameters(blk);
  return std::make_shared<LuaSpatialMaterialFunction>(params);
}

} // namespace

RegisterLuaFunctionAsIs(chiAdjointSolverAddResponseFunction);

int
chiAdjointSolverAddResponseFunction(lua_State* L)
{
  const std::string fname = "chiAdjointSolverAddResponseFunction";
  const int num_args = lua_gettop(L);
  if (num_args < 3) LuaPostArgAmountError(fname, 3, num_args);

  LuaCheckNilValue(fname, L, 1);
  LuaCheckNilValue(fname, L, 2);
  LuaCheckNilValue(fname, L, 3);

  LuaCheckIntegerValue(fname, L, 1);
  LuaCheckStringValue(fname, L, 2);
  LuaCheckNumberValue(fname, L, 3);

  const int solver_handle = lua_tointeger(L, 1);
  const std::string qoi_name = lua_tostring(L, 2);
  const int logvol_handle = lua_tointeger(L, 3);

  std::string lua_function;
  if (num_args == 4)
  {
    LuaCheckNilValue(fname, L, 4);

    lua_function = lua_tostring(L, 4);
  }

  auto response_function = CreateResponseFunction(lua_function);
  opensn::Chi::function_stack.push_back(response_function);

  auto& solver = opensn::Chi::GetStackItem<opensn::lbs::DiscreteOrdinatesAdjointSolver>(
    opensn::Chi::object_stack, solver_handle, fname);

  auto p_logical_volume = std::dynamic_pointer_cast<LogicalVolume>(
    opensn::Chi::GetStackItemPtr(opensn::Chi::object_stack, logvol_handle, fname));

  size_t qoi_index = solver.AddResponseFunction(qoi_name, p_logical_volume, response_function);
  lua_pushinteger(L, static_cast<lua_Integer>(qoi_index));

  return 1;
}

} // namespace opensnlua::lbs
