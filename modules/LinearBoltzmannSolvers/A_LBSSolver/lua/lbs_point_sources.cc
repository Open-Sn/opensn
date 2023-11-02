#include "A_LBSSolver/lbs_solver.h"

#include "chi_runtime.h"
#include "chi_log.h"

namespace lbs::common_lua_utils
{

//###################################################################
/**Adds a point source to an LBS solver.
\param SolverIndex int Handle to the solver.
\param Location_x double X-location.
\param Location_y double Y-location.
\param Location_z double Z-location.
\param Strength table Source strength as a multigroup vector.

 \ingroup LBSLuaFunctions
 */
int chiLBSAddPointSource(lua_State *L)
{
  const std::string fname = "chiLBSAddPointSource";
  const int num_args = lua_gettop(L);
  if (num_args != 5)
    LuaPostArgAmountError(fname, 5, num_args);

  LuaCheckNilValue(fname, L, 1);
  LuaCheckNilValue(fname, L, 2);
  LuaCheckNilValue(fname, L, 3);
  LuaCheckNilValue(fname, L, 4);
  LuaCheckNilValue(fname, L, 5);

  //============================================= Get pointer to solver
  const int solver_handle = lua_tonumber(L,1);
  auto& lbs_solver =
    Chi::GetStackItem<lbs::LBSSolver>(Chi::object_stack,
                                                       solver_handle,
                                                       fname);

  //============================================= Get other arguments
  const double x = lua_tonumber(L, 2);
  const double y = lua_tonumber(L, 3);
  const double z = lua_tonumber(L, 4);

  const chi_mesh::Vector3 location(x,y,z);

  LuaCheckTableValue(fname, L, 5);

  std::vector<double> groupwise_strength;
  LuaPopulateVectorFrom1DArray(fname, L, 5, groupwise_strength);

  lbs_solver.AddPointSource(lbs::PointSource(location, groupwise_strength));

  Chi::log.Log() << "LBS: Added point source at "
                << location.PrintStr();

  return 0;
}

//###################################################################
/**Clears all the point sources from the solver. This is mostly
 * useful for adjoint response calculations.
\param SolverIndex int Handle to the solver.

 \ingroup LBSLuaFunctions
 */
int chiLBSClearPointSources(lua_State *L)
{
  const std::string fname = "chiLBSClearPointSources";
  const int num_args = lua_gettop(L);
  if (num_args != 1)
    LuaPostArgAmountError(fname, 1, num_args);

  LuaCheckNilValue(fname, L, 1);

  //============================================= Get pointer to solver
  const int solver_handle = lua_tonumber(L,1);
  auto& lbs_solver =
    Chi::GetStackItem<lbs::LBSSolver>(Chi::object_stack,
                                                       solver_handle,
                                                       fname);

  lbs_solver.ClearPointSources();

  Chi::log.Log() << "LBS: Cleared all point sources";

  return 0;
}

//###################################################################
/**Initializes the point sources. This is mostly
 * useful for adjoint response calculations.
\param SolverIndex int Handle to the solver.

 \ingroup LBSLuaFunctions
 */
int chiLBSInitializePointSources(lua_State *L)
{
  const std::string fname = "chiLBSInitializePointSources";
  const int num_args = lua_gettop(L);
  if (num_args != 1)
    LuaPostArgAmountError(fname, 1, num_args);

  LuaCheckNilValue(fname, L, 1);

  //============================================= Get pointer to solver
  const int solver_handle = lua_tonumber(L,1);
  auto& lbs_solver =
    Chi::GetStackItem<lbs::LBSSolver>(Chi::object_stack,
                                                       solver_handle,
                                                       fname);

  lbs_solver.InitializePointSources();

  Chi::log.Log() << "LBS: Initializing point sources.";

  return 0;
}

}//namespace lbs::common_lua_utils