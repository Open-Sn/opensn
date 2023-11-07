#include "modules/LinearBoltzmannSolvers/A_LBSSolver/lbs_solver.h"

#include "framework/runtime.h"
#include "framework/logging/log.h"

namespace lbs::common_lua_utils
{

int
chiLBSWriteFluxMoments(lua_State* L)
{
  const std::string fname = "chiLBSWriteFluxMoments";
  // Get arguments
  const int num_args = lua_gettop(L);
  if (num_args != 2) LuaPostArgAmountError(fname, 2, num_args);

  LuaCheckNilValue(fname, L, 1);
  LuaCheckNilValue(fname, L, 2);

  const int solver_handle = lua_tonumber(L, 1);
  const std::string file_base = lua_tostring(L, 2);

  // Get pointer to solver
  auto& lbs_solver = Chi::GetStackItem<lbs::LBSSolver>(Chi::object_stack, solver_handle, fname);

  lbs_solver.WriteFluxMoments(file_base, lbs_solver.PhiOldLocal());

  return 0;
}

int
chiLBSCreateAndWriteSourceMoments(lua_State* L)
{
  const std::string fname = "chiLBSCreateAndWriteSourceMoments";
  // Get arguments
  const int num_args = lua_gettop(L);
  if (num_args != 2) LuaPostArgAmountError(fname, 2, num_args);

  LuaCheckNilValue(fname, L, 1);
  LuaCheckNilValue(fname, L, 2);

  const int solver_handle = lua_tonumber(L, 1);
  const std::string file_base = lua_tostring(L, 2);

  // Get pointer to solver
  auto& lbs_solver = Chi::GetStackItem<lbs::LBSSolver>(Chi::object_stack, solver_handle, fname);

  auto source_moments = lbs_solver.MakeSourceMomentsFromPhi();
  lbs_solver.WriteFluxMoments(file_base, source_moments);

  return 0;
}

int
chiLBSReadFluxMomentsAndMakeSourceMoments(lua_State* L)
{
  const std::string fname = "chiLBSReadFluxMomentsAndMakeSourceMoments";
  // Get arguments
  const int num_args = lua_gettop(L);
  if ((num_args != 2) and (num_args != 3)) LuaPostArgAmountError(fname, 2, num_args);

  LuaCheckNilValue(fname, L, 1);
  LuaCheckNilValue(fname, L, 2);

  const int solver_handle = lua_tonumber(L, 1);
  const std::string file_base = lua_tostring(L, 2);

  bool single_file_flag = false;
  if (num_args == 3)
  {
    LuaCheckBoolValue(fname, L, 3);
    single_file_flag = lua_toboolean(L, 3);
  }

  // Get pointer to solver
  auto& lbs_solver = Chi::GetStackItem<lbs::LBSSolver>(Chi::object_stack, solver_handle, fname);

  lbs_solver.ReadFluxMoments(file_base, lbs_solver.ExtSrcMomentsLocal(), single_file_flag);

  Chi::log.Log() << "Making source moments from flux file.";
  auto temp_phi = lbs_solver.PhiOldLocal();
  lbs_solver.PhiOldLocal() = lbs_solver.ExtSrcMomentsLocal();
  lbs_solver.ExtSrcMomentsLocal() = lbs_solver.MakeSourceMomentsFromPhi();
  lbs_solver.PhiOldLocal() = temp_phi;

  return 0;
}

int
chiLBSReadSourceMoments(lua_State* L)
{
  const std::string fname = "chiLBSReadSourceMoments";
  // Get arguments
  const int num_args = lua_gettop(L);
  if ((num_args != 2) and (num_args != 3)) LuaPostArgAmountError(fname, 2, num_args);

  LuaCheckNilValue(fname, L, 1);
  LuaCheckNilValue(fname, L, 2);

  const int solver_handle = lua_tonumber(L, 1);
  const std::string file_base = lua_tostring(L, 2);

  bool single_file_flag = false;
  if (num_args == 3)
  {
    LuaCheckBoolValue(fname, L, 3);
    single_file_flag = lua_toboolean(L, 3);
  }

  // Get pointer to solver
  auto& lbs_solver = Chi::GetStackItem<lbs::LBSSolver>(Chi::object_stack, solver_handle, fname);

  lbs_solver.ReadFluxMoments(file_base, lbs_solver.ExtSrcMomentsLocal(), single_file_flag);

  return 0;
}

int
chiLBSReadFluxMoments(lua_State* L)
{
  const std::string fname = "chiLBSReadFluxMoments";
  // Get arguments
  const int num_args = lua_gettop(L);
  if ((num_args != 2) and (num_args != 3)) LuaPostArgAmountError(fname, 2, num_args);

  LuaCheckNilValue(fname, L, 1);
  LuaCheckNilValue(fname, L, 2);

  const int solver_handle = lua_tonumber(L, 1);
  const std::string file_base = lua_tostring(L, 2);

  bool single_file_flag = false;
  if (num_args == 3)
  {
    LuaCheckBoolValue(fname, L, 3);
    single_file_flag = lua_toboolean(L, 3);
  }

  // Get pointer to solver
  auto& lbs_solver = Chi::GetStackItem<lbs::LBSSolver>(Chi::object_stack, solver_handle, fname);

  lbs_solver.ReadFluxMoments(file_base, lbs_solver.PhiOldLocal(), single_file_flag);

  return 0;
}

} // namespace lbs::common_lua_utils
