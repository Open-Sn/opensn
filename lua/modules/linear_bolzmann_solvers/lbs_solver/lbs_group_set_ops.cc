#include "modules/linear_boltzmann_solvers/a_lbs_solver/lbs_solver.h"
#include "framework/math/quadratures/angular_product_quadrature.h"
#include "modules/linear_boltzmann_solvers/a_lbs_solver/groupset/lbs_group.h"
#include "modules/linear_boltzmann_solvers/a_lbs_solver/groupset/lbs_groupset.h"
#include "framework/runtime.h"
#include "framework/logging/log.h"
#include "framework/lua.h"

using namespace opensn;

namespace opensnlua::lbs
{

int
LBSCreateGroupset(lua_State* L)
{
  const std::string fname = "LBSCreateGroupset";
  // Get pointer to solver
  const int solver_handle = lua_tonumber(L, 1);
  auto& lbs_solver =
    opensn::GetStackItem<opensn::lbs::LBSSolver>(opensn::object_stack, solver_handle, fname);

  // Create groupset
  lbs_solver.AddGroupset();

  lua_pushinteger(L, lbs_solver.Groupsets().back().id_);
  return 1;
}

int
LBSCreateGroup(lua_State* L)
{
  const std::string fname = "LBSCreateGroup";
  const int num_args = lua_gettop(L);
  if (num_args < 1) LuaPostArgAmountError(fname, 1, num_args);

  // Get solver
  LuaCheckNumberValue(fname, L, 1);
  const int solver_handle = lua_tointeger(L, 1);

  auto& lbs_solver =
    opensn::GetStackItem<opensn::lbs::LBSSolver>(opensn::object_stack, solver_handle, fname);

  // Set Id
  int group_id = -1;
  if (num_args == 2)
  {
    LuaCheckNumberValue(fname, L, 2);
    group_id = lua_tointeger(L, 2);
  }

  // Create groupset
  lbs_solver.AddGroup(group_id);

  lua_pushinteger(L, lbs_solver.Groups().back().id_);
  return 1;
}

int
LBSGroupsetAddGroups(lua_State* L)
{
  const std::string fname = "LBSGroupsetAddGroups";
  // Get arguments
  const int num_args = lua_gettop(L);
  if (num_args != 4) LuaPostArgAmountError(fname, 4, num_args);
  int solver_handle = lua_tonumber(L, 1);
  int grpset_index = lua_tonumber(L, 2);
  int from = lua_tonumber(L, 3);
  int to = lua_tonumber(L, 4);

  // Get pointer to solver
  auto& lbs_solver =
    opensn::GetStackItem<opensn::lbs::LBSSolver>(opensn::object_stack, solver_handle, fname);

  // Obtain pointer to groupset
  opensn::lbs::LBSGroupset* groupset = nullptr;
  try
  {
    groupset = &lbs_solver.Groupsets().at(grpset_index);
  }
  catch (const std::out_of_range& o)
  {
    opensn::log.LogAllError() << "LBSGroupsetAddGroups: Invalid handle to groupset\n";
    opensn::Exit(EXIT_FAILURE);
  }
  if (groupset == nullptr) throw std::runtime_error("LBSGroupsetAddGroups: Bad trouble.");

  // Add the groups
  if (to < from)
  {
    opensn::log.LogAllError() << "No groups added to groupset in LBSGroupsetAddGroups. "
                                 "This is triggered when groups are added with the \"to\" "
                                 "field being less than the \"from\" field.";
    opensn::Exit(EXIT_FAILURE);
  }

  for (unsigned k = from; k <= to; k++)
  {
    opensn::lbs::LBSGroup const* group = nullptr;
    // Check valid group
    try
    {
      group = &lbs_solver.Groups().at(k);
    }
    catch (const std::out_of_range& o)
    {
      opensn::log.LogAllError() << "LBSGroupsetAddGroups: Invalid group added to groupset\n";
      opensn::Exit(EXIT_FAILURE);
    }
    if (group == nullptr) throw std::runtime_error("LBSGroupsetAddGroups: Bad trouble.");

    groupset->groups_.push_back(*group);
  }
  return 0;
}

int
LBSGroupsetSetQuadrature(lua_State* L)
{
  const std::string fname = "LBSGroupsetSetQuadrature";
  // Get arguments
  const int num_args = lua_gettop(L);
  if (num_args != 3) LuaPostArgAmountError(fname, 3, num_args);

  LuaCheckNilValue(fname, L, 1);
  LuaCheckNilValue(fname, L, 2);
  LuaCheckNilValue(fname, L, 3);

  int solver_handle = lua_tonumber(L, 1);
  int grpset_index = lua_tonumber(L, 2);
  int prquad_index = lua_tonumber(L, 3);

  // Get pointer to solver
  auto& lbs_solver =
    opensn::GetStackItem<opensn::lbs::LBSSolver>(opensn::object_stack, solver_handle, fname);

  // Obtain pointer to groupset
  opensn::lbs::LBSGroupset* groupset = nullptr;
  try
  {
    groupset = &lbs_solver.Groupsets().at(grpset_index);
  }
  catch (const std::out_of_range& o)
  {
    opensn::log.LogAllError() << "Invalid handle to groupset"
                                 "in LBSGroupsetSetQuadrature.";
    opensn::Exit(EXIT_FAILURE);
  }
  if (groupset == nullptr) throw std::logic_error("LBSGroupsetSetQuadrature: Bad trouble");

  // Obtain pointer to quadrature
  std::shared_ptr<AngularQuadrature> ang_quad;
  try
  {
    ang_quad = opensn::GetStackItemPtr(opensn::angular_quadrature_stack, prquad_index);
  }
  catch (const std::out_of_range& o)
  {
    opensn::log.LogAllError() << "Invalid handle to Product Quadrature"
                                 "in LBSGroupsetSetQuadrature. Handle provided: "
                              << prquad_index;
    opensn::Exit(EXIT_FAILURE);
  }

  groupset->quadrature_ = ang_quad;

  if (ang_quad->type_ == AngularQuadratureType::ProductQuadrature)
  {
    auto prodquad = std::static_pointer_cast<ProductQuadrature>(ang_quad);

    opensn::log.Log() << "Groupset " << grpset_index << " quadrature set to quadrature with "
                      << prodquad->azimu_ang_.size() << " azimuthal angles and "
                      << prodquad->polar_ang_.size() << " polar angles. ";
  }
  else if (ang_quad->type_ == AngularQuadratureType::Arbitrary)
  {
    opensn::log.Log() << "Groupset " << grpset_index << " quadrature set to quadrature with "
                      << ang_quad->abscissae_.size() << " number of angles.";
  }
  else
    opensn::log.Log() << "Groupset " << grpset_index << " quadrature set unknown quadrature type";

  return 0;
}

int
LBSGroupsetSetAngleAggregationType(lua_State* L)
{
  const std::string fname = "LBSGroupsetSetAngleAggregationType";
  // Get arguments
  const int num_args = lua_gettop(L);
  if (num_args != 3) LuaPostArgAmountError(fname, 3, num_args);

  LuaCheckNilValue(fname, L, 1);
  LuaCheckNilValue(fname, L, 2);
  LuaCheckNilValue(fname, L, 3);
  int solver_handle = lua_tonumber(L, 1);
  int grpset_index = lua_tonumber(L, 2);
  int agg_type = lua_tonumber(L, 3);

  // Get pointer to solver
  auto& lbs_solver =
    opensn::GetStackItem<opensn::lbs::LBSSolver>(opensn::object_stack, solver_handle, fname);

  // Obtain pointer to groupset
  opensn::lbs::LBSGroupset* groupset = nullptr;
  try
  {
    groupset = &lbs_solver.Groupsets().at(grpset_index);
  }
  catch (const std::out_of_range& o)
  {
    opensn::log.LogAllError() << "Invalid handle to groupset "
                              << "in call to LBSGroupsetSetAngleAggregationType";
    opensn::Exit(EXIT_FAILURE);
  }

  // Setting aggregation type
  if (agg_type == (int)opensn::lbs::AngleAggregationType::SINGLE)
    groupset->angleagg_method_ = opensn::lbs::AngleAggregationType::SINGLE;
  else if (agg_type == (int)opensn::lbs::AngleAggregationType::POLAR)
    groupset->angleagg_method_ = opensn::lbs::AngleAggregationType::POLAR;
  else if (agg_type == (int)opensn::lbs::AngleAggregationType::AZIMUTHAL)
    groupset->angleagg_method_ = opensn::lbs::AngleAggregationType::AZIMUTHAL;
  else
  {
    opensn::log.LogAllError() << "Invalid aggregation type to groupset " << grpset_index
                              << " in call to LBSGroupsetSetAngleAggregationType";
    opensn::Exit(EXIT_FAILURE);
  }

  opensn::log.Log() << "Groupset " << grpset_index << " Angle aggregation set to " << agg_type;

  return 0;
}

int
LBSGroupsetSetAngleAggDiv(lua_State* L)
{
  const std::string fname = "LBSGroupsetSetAngleAggDiv";
  // Get arguments
  const int num_args = lua_gettop(L);
  if (num_args != 3) LuaPostArgAmountError(fname, 3, num_args);

  LuaCheckNilValue(fname, L, 1);
  LuaCheckNilValue(fname, L, 2);
  LuaCheckNilValue(fname, L, 3);
  int solver_handle = lua_tonumber(L, 1);
  int grpset_index = lua_tonumber(L, 2);
  int num_div = lua_tonumber(L, 3);

  // Get pointer to solver
  auto& lbs_solver =
    opensn::GetStackItem<opensn::lbs::LBSSolver>(opensn::object_stack, solver_handle, fname);

  // Obtain pointer to groupset
  opensn::lbs::LBSGroupset* groupset = nullptr;
  try
  {
    groupset = &lbs_solver.Groupsets().at(grpset_index);
  }
  catch (const std::out_of_range& o)
  {
    opensn::log.LogAllError() << "Invalid handle to groupset "
                              << "in call to LBSGroupsetSetAngleAggDiv";
    opensn::Exit(EXIT_FAILURE);
  }

  // Bounds checking
  if (num_div <= 0)
  {
    opensn::log.LogAllError() << "Invalid number of divisions "
                              << "in call to LBSGroupsetSetAngleAggDiv. Must be >= 1.";
  }

  groupset->master_num_ang_subsets_ = num_div;

  opensn::log.Log() << "Groupset " << grpset_index << " angle aggregation divisions "
                    << "set to " << num_div;

  return 0;
}

int
LBSGroupsetSetGroupSubsets(lua_State* L)
{
  const std::string fname = "LBSGroupsetSetGroupSubsets";
  // Get arguments
  const int num_args = lua_gettop(L);
  if (num_args != 3) LuaPostArgAmountError(fname, 3, num_args);

  LuaCheckNilValue(fname, L, 1);
  LuaCheckNilValue(fname, L, 2);
  LuaCheckNilValue(fname, L, 3);
  int solver_handle = lua_tonumber(L, 1);
  int grpset_index = lua_tonumber(L, 2);
  int num_div = lua_tonumber(L, 3);

  // Get pointer to solver
  auto& lbs_solver =
    opensn::GetStackItem<opensn::lbs::LBSSolver>(opensn::object_stack, solver_handle, fname);

  // Obtain pointer to groupset
  opensn::lbs::LBSGroupset* groupset = nullptr;
  try
  {
    groupset = &lbs_solver.Groupsets().at(grpset_index);
  }
  catch (const std::out_of_range& o)
  {
    opensn::log.LogAllError() << "Invalid handle to groupset "
                              << "in call to LBSGroupsetSetGroupSubsets";
    opensn::Exit(EXIT_FAILURE);
  }

  // Bounds checking
  if (num_div <= 0)
  {
    opensn::log.LogAllError() << "Invalid number of subsets "
                              << "in call to LBSGroupsetSetGroupSubsets. Must be >= 1.";
  }

  groupset->master_num_grp_subsets_ = num_div;

  opensn::log.Log() << "Groupset " << grpset_index << " subset divisions "
                    << "set to " << num_div;

  return 0;
}

int
LBSGroupsetSetIterativeMethod(lua_State* L)
{
  const std::string fname = "LBSGroupsetSetIterativeMethod";
  // Get arguments
  const int num_args = lua_gettop(L);
  if (num_args != 3) LuaPostArgAmountError(fname, 3, num_args);

  LuaCheckNilValue(fname, L, 1);
  LuaCheckNilValue(fname, L, 2);
  LuaCheckNilValue(fname, L, 3);
  int solver_handle = lua_tonumber(L, 1);
  int grpset_index = lua_tonumber(L, 2);
  int iter_method = lua_tonumber(L, 3);

  // Get pointer to solver
  auto& lbs_solver =
    opensn::GetStackItem<opensn::lbs::LBSSolver>(opensn::object_stack, solver_handle, fname);

  // Obtain pointer to groupset
  opensn::lbs::LBSGroupset* groupset = nullptr;
  try
  {
    groupset = &lbs_solver.Groupsets().at(grpset_index);
  }
  catch (const std::out_of_range& o)
  {
    opensn::log.LogAllError() << "Invalid handle to groupset "
                              << "in call to LBSGroupsetSetGroupSubsets";
    opensn::Exit(EXIT_FAILURE);
  }

  {
    using opensn::lbs::IterativeMethod;
    if (iter_method == static_cast<int>(IterativeMethod::CLASSICRICHARDSON))
    {
      groupset->iterative_method_ = IterativeMethod::CLASSICRICHARDSON;
    }
    else if (iter_method == static_cast<int>(IterativeMethod::CLASSICRICHARDSON_CYCLES))
    {
      groupset->allow_cycles_ = true;
      groupset->iterative_method_ = IterativeMethod::CLASSICRICHARDSON;
    }
    else if (iter_method == static_cast<int>(IterativeMethod::GMRES))
    {
      throw std::invalid_argument(fname + "Deprecated iterative method GMRES, "
                                          "use KRYLOV_GMRES.");
    }
    else if (iter_method == static_cast<int>(IterativeMethod::GMRES_CYCLES))
    {
      throw std::invalid_argument(fname + "Deprecated iterative method GMRES_CYCLES, "
                                          "use KRYLOV_GMRES_CYCLES.");
    }
    else if (iter_method == static_cast<int>(IterativeMethod::KRYLOV_RICHARDSON))
    {
      groupset->iterative_method_ = IterativeMethod::KRYLOV_RICHARDSON;
    }
    else if (iter_method == static_cast<int>(IterativeMethod::KRYLOV_RICHARDSON_CYCLES))
    {
      groupset->allow_cycles_ = true;
      groupset->iterative_method_ = IterativeMethod::KRYLOV_RICHARDSON;
    }
    else if (iter_method == static_cast<int>(IterativeMethod::KRYLOV_GMRES))
    {
      groupset->iterative_method_ = IterativeMethod::KRYLOV_GMRES;
    }
    else if (iter_method == static_cast<int>(IterativeMethod::KRYLOV_GMRES_CYCLES))
    {
      groupset->allow_cycles_ = true;
      groupset->iterative_method_ = IterativeMethod::KRYLOV_GMRES;
    }
    else if (iter_method == static_cast<int>(IterativeMethod::KRYLOV_BICGSTAB))
    {
      groupset->iterative_method_ = IterativeMethod::KRYLOV_BICGSTAB;
    }
    else if (iter_method == static_cast<int>(IterativeMethod::KRYLOV_BICGSTAB_CYCLES))
    {
      groupset->allow_cycles_ = true;
      groupset->iterative_method_ = IterativeMethod::KRYLOV_BICGSTAB;
    }
    else
    {
      opensn::log.LogAllError() << "Unsupported iterative method specified in call to "
                                << "LBSGroupsetSetIterativeMethod.";
      opensn::Exit(EXIT_FAILURE);
    }
  }

  return 0;
}

int
LBSGroupsetSetResidualTolerance(lua_State* L)
{
  const std::string fname = "LBSGroupsetSetResidualTolerance";
  // Get arguments
  const int num_args = lua_gettop(L);
  if (num_args != 3) LuaPostArgAmountError(fname, 3, num_args);

  LuaCheckNilValue(fname, L, 1);
  LuaCheckNilValue(fname, L, 2);
  LuaCheckNilValue(fname, L, 3);
  int solver_handle = lua_tonumber(L, 1);
  int grpset_index = lua_tonumber(L, 2);
  double resid_tol = lua_tonumber(L, 3);

  // Get pointer to solver
  auto& lbs_solver =
    opensn::GetStackItem<opensn::lbs::LBSSolver>(opensn::object_stack, solver_handle, fname);

  // Obtain pointer to groupset
  opensn::lbs::LBSGroupset* groupset = nullptr;
  try
  {
    groupset = &lbs_solver.Groupsets().at(grpset_index);
  }
  catch (const std::out_of_range& o)
  {
    opensn::log.LogAllError() << "Invalid handle to groupset "
                              << "in call to LBSGroupsetSetGroupSubsets";
    opensn::Exit(EXIT_FAILURE);
  }

  // Bounds checking
  if (resid_tol < 0)
  {
    opensn::log.LogAllError() << "Invalid residual tolerance specified. Must be greater >= 0.0";
    opensn::Exit(EXIT_FAILURE);
  }

  groupset->residual_tolerance_ = resid_tol;

  char buff[100];
  snprintf(buff, 100, "%.4e", resid_tol);

  opensn::log.Log() << "Groupset " << grpset_index << " residual tolerance "
                    << "set to " << buff;

  return 0;
}

int
LBSGroupsetSetMaxIterations(lua_State* L)
{
  const std::string fname = "LBSGroupsetSetMaxIterations";
  // Get arguments
  const int num_args = lua_gettop(L);
  if (num_args != 3) LuaPostArgAmountError(fname, 3, num_args);

  LuaCheckNilValue(fname, L, 1);
  LuaCheckNilValue(fname, L, 2);
  LuaCheckNilValue(fname, L, 3);
  int solver_handle = lua_tonumber(L, 1);
  int grpset_index = lua_tonumber(L, 2);
  int num_iter = lua_tonumber(L, 3);

  // Get pointer to solver
  auto& lbs_solver =
    opensn::GetStackItem<opensn::lbs::LBSSolver>(opensn::object_stack, solver_handle, fname);

  // Obtain pointer to groupset
  opensn::lbs::LBSGroupset* groupset = nullptr;
  try
  {
    groupset = &lbs_solver.Groupsets().at(grpset_index);
  }
  catch (const std::out_of_range& o)
  {
    opensn::log.LogAllError() << "Invalid handle to groupset "
                              << "in call to LBSGroupsetSetMaxIterations";
    opensn::Exit(EXIT_FAILURE);
  }

  // Bounds checking
  if (num_iter < 0)
  {
    opensn::log.LogAllError() << "Invalid number of iterations "
                              << "in call to LBSGroupsetSetMaxIterations. Must be >= 0.";
  }

  groupset->max_iterations_ = num_iter;

  opensn::log.Log() << "Groupset " << grpset_index << " max # iterations "
                    << "set to " << num_iter;

  return 0;
}

int
LBSGroupsetSetGMRESRestartIntvl(lua_State* L)
{
  const std::string fname = "LBSGroupsetSetGMRESRestartIntvl";
  // Get arguments
  const int num_args = lua_gettop(L);
  if (num_args != 3) LuaPostArgAmountError(fname, 3, num_args);

  LuaCheckNilValue(fname, L, 1);
  LuaCheckNilValue(fname, L, 2);
  LuaCheckNilValue(fname, L, 3);
  int solver_handle = lua_tonumber(L, 1);
  int grpset_index = lua_tonumber(L, 2);
  int restart_intvl = lua_tonumber(L, 3);

  // Get pointer to solver
  auto& lbs_solver =
    opensn::GetStackItem<opensn::lbs::LBSSolver>(opensn::object_stack, solver_handle, fname);

  // Obtain pointer to groupset
  opensn::lbs::LBSGroupset* groupset = nullptr;
  try
  {
    groupset = &lbs_solver.Groupsets().at(grpset_index);
  }
  catch (const std::out_of_range& o)
  {
    opensn::log.LogAllError() << "Invalid handle to groupset "
                              << "in call to LBSGroupsetSetGMRESRestartIntvl";
    opensn::Exit(EXIT_FAILURE);
  }

  // Bounds checking
  if (restart_intvl < 2)
  {
    opensn::log.LogAllError() << "Invalid GMRES restart interval specified "
                              << "in call to LBSGroupsetSetGMRESRestartIntvl. Must be >= 3.";
  }

  groupset->gmres_restart_intvl_ = restart_intvl;

  opensn::log.Log() << "Groupset " << grpset_index << " GMRES restart interval set to "
                    << "set to " << restart_intvl;

  return 0;
}

int
LBSGroupsetSetEnableSweepLog(lua_State* L)
{
  const std::string fname = "LBSGroupsetSetEnableSweepLog";
  // Get arguments
  const int num_args = lua_gettop(L);
  if (num_args != 3) LuaPostArgAmountError(fname, 3, num_args);

  LuaCheckNilValue(fname, L, 1);
  LuaCheckNilValue(fname, L, 2);
  LuaCheckNilValue(fname, L, 3);
  int solver_handle = lua_tonumber(L, 1);
  int grpset_index = lua_tonumber(L, 2);
  bool log_flag = lua_toboolean(L, 3);

  // Get pointer to solver
  auto& lbs_solver =
    opensn::GetStackItem<opensn::lbs::LBSSolver>(opensn::object_stack, solver_handle, fname);

  // Obtain pointer to groupset
  opensn::lbs::LBSGroupset* groupset = nullptr;
  try
  {
    groupset = &lbs_solver.Groupsets().at(grpset_index);
  }
  catch (const std::out_of_range& o)
  {
    opensn::log.LogAllError() << "Invalid handle to groupset "
                              << "in call to LBSGroupsetSetEnableSweepLog";
    opensn::Exit(EXIT_FAILURE);
  }

  groupset->log_sweep_events_ = log_flag;

  opensn::log.Log() << "Groupset " << grpset_index << " flag for writing sweep log "
                    << "set to " << log_flag;

  return 0;
}

int
LBSGroupsetSetWGDSA(lua_State* L)
{
  const std::string fname = "LBSGroupsetSetWGDSA";
  // Get arguments
  const int num_args = lua_gettop(L);
  if (num_args < 4) LuaPostArgAmountError(fname, 4, num_args);

  LuaCheckNilValue(fname, L, 1);
  LuaCheckNilValue(fname, L, 2);
  LuaCheckNilValue(fname, L, 3);
  LuaCheckNilValue(fname, L, 4);
  int solver_handle = lua_tonumber(L, 1);
  int grpset_index = lua_tonumber(L, 2);
  int max_iters = lua_tonumber(L, 3);
  double resid_tol = lua_tonumber(L, 4);
  bool verbose = false;
  const char* petsc_string = "";

  if (num_args >= 5) verbose = lua_toboolean(L, 5);

  if (num_args == 6) petsc_string = lua_tostring(L, 6);

  // Get pointer to solver
  auto& lbs_solver =
    opensn::GetStackItem<opensn::lbs::LBSSolver>(opensn::object_stack, solver_handle, fname);

  // Obtain pointer to groupset
  opensn::lbs::LBSGroupset* groupset = nullptr;
  try
  {
    groupset = &lbs_solver.Groupsets().at(grpset_index);
  }
  catch (const std::out_of_range& o)
  {
    opensn::log.LogAllError() << "Invalid handle to groupset "
                              << "in call to LBSGroupsetSetWGDSA";
    opensn::Exit(EXIT_FAILURE);
  }

  groupset->apply_wgdsa_ = true;
  groupset->wgdsa_max_iters_ = max_iters;
  groupset->wgdsa_tol_ = resid_tol;
  groupset->wgdsa_verbose_ = verbose;
  groupset->wgdsa_string_ = std::string(petsc_string);

  opensn::log.Log() << "Groupset " << grpset_index << " set to apply WGDSA with " << max_iters
                    << " maximum iterations and a tolerance of " << resid_tol
                    << ". PETSc-string: " << std::string(petsc_string);

  return 0;
}

int
LBSGroupsetSetTGDSA(lua_State* L)
{
  const std::string fname = "LBSGroupsetSetTGDSA";
  // Get arguments
  const int num_args = lua_gettop(L);
  if (num_args < 4) LuaPostArgAmountError(fname, 4, num_args);

  LuaCheckNilValue(fname, L, 1);
  LuaCheckNilValue(fname, L, 2);
  LuaCheckNilValue(fname, L, 3);
  LuaCheckNilValue(fname, L, 4);
  int solver_handle = lua_tonumber(L, 1);
  int grpset_index = lua_tonumber(L, 2);
  int max_iters = lua_tonumber(L, 3);
  double resid_tol = lua_tonumber(L, 4);
  bool verbose = false;
  const char* petsc_string = "";

  if (num_args >= 5) verbose = lua_toboolean(L, 5);

  if (num_args == 6) petsc_string = lua_tostring(L, 6);

  // Get pointer to solver
  auto& lbs_solver =
    opensn::GetStackItem<opensn::lbs::LBSSolver>(opensn::object_stack, solver_handle, fname);

  // Obtain pointer to groupset
  opensn::lbs::LBSGroupset* groupset = nullptr;
  try
  {
    groupset = &lbs_solver.Groupsets().at(grpset_index);
  }
  catch (const std::out_of_range& o)
  {
    opensn::log.LogAllError() << "Invalid handle to groupset "
                              << "in call to LBSGroupsetSetTGDSA";
    opensn::Exit(EXIT_FAILURE);
  }

  groupset->apply_tgdsa_ = true;
  groupset->tgdsa_max_iters_ = max_iters;
  groupset->tgdsa_tol_ = resid_tol;
  groupset->tgdsa_verbose_ = verbose;
  groupset->tgdsa_string_ = std::string(petsc_string);

  opensn::log.Log() << "Groupset " << grpset_index << " set to apply TGDSA with " << max_iters
                    << " maximum iterations and a tolerance of " << resid_tol
                    << ". PETSc-string: " << std::string(petsc_string);

  return 0;
}

} // namespace opensnlua::lbs
