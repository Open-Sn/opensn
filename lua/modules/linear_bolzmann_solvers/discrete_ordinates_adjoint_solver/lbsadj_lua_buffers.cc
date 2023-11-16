#include "lbsadj_lua_utils.h"

#include "modules/linear_boltzmann_solvers/c_discrete_ordinates_adjoint_solver/lbs_adj_solver.h"

#include "framework/runtime.h"

#include "framework/console/console.h"

namespace opensnlua::lbs
{

RegisterLuaFunctionAsIs(chiAdjointSolverReadFluxMomentsToBuffer);
RegisterLuaFunctionAsIs(chiAdjointSolverApplyFluxMomentBuffer);

int
chiAdjointSolverReadFluxMomentsToBuffer(lua_State* L)
{
  const std::string fname = __FUNCTION__;
  const int num_args = lua_gettop(L);
  if (num_args != 2) LuaPostArgAmountError(fname, 2, num_args);

  LuaCheckNilValue(fname, L, 1);
  LuaCheckNilValue(fname, L, 2);

  const int solver_handle = lua_tointeger(L, 1);

  auto& solver = opensn::Chi::GetStackItem<opensn::lbs::DiscreteOrdinatesAdjointSolver>(
    opensn::Chi::object_stack, solver_handle, fname);

  const std::string file_basename = lua_tostring(L, 2);

  std::vector<double> moments;
  solver.ReadFluxMoments(file_basename, moments);

  solver.m_moment_buffers_.push_back(std::move(moments));

  const size_t handle = solver.m_moment_buffers_.size() - 1;
  lua_pushinteger(L, static_cast<lua_Integer>(handle));

  return 1;
}

int
chiAdjointSolverApplyFluxMomentBuffer(lua_State* L)
{
  const std::string fname = __FUNCTION__;
  const int num_args = lua_gettop(L);
  if (num_args != 2) LuaPostArgAmountError(fname, 2, num_args);

  LuaCheckNilValue(fname, L, 1);
  LuaCheckNilValue(fname, L, 2);

  const int solver_handle = lua_tointeger(L, 1);

  auto& solver = opensn::Chi::GetStackItem<opensn::lbs::DiscreteOrdinatesAdjointSolver>(
    opensn::Chi::object_stack, solver_handle, fname);

  const int buffer_handle = lua_tointeger(L, 2);

  if (buffer_handle < 0 or buffer_handle >= solver.m_moment_buffers_.size())
    throw std::invalid_argument(fname + ": Invalid buffer handle.");

  solver.PhiOldLocal() = solver.m_moment_buffers_[buffer_handle];

  return 0;
}

} // namespace opensnlua::lbs
