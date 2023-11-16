#include "modules/linear_boltzmann_solvers/a_lbs_solver/lbs_solver.h"

#define scpcode static_cast<PropertyCode>
#define scint static_cast<int>

#include "framework/runtime.h"
#include "framework/logging/log.h"

using namespace opensn;

namespace opensnlua::lbs
{

enum class PropertyCode : int
{
  DISCRETIZATION_METHOD = 1,
  PWLD = 3,
  BOUNDARY_CONDITION = 3,
  XMAX = 31,
  XMIN = 32,
  YMAX = 33,
  YMIN = 34,
  ZMAX = 35,
  ZMIN = 36,
  SCATTERING_ORDER = 4,
  SWEEP_EAGER_LIMIT = 5,
  READ_RESTART_DATA = 6,
  WRITE_RESTART_DATA = 7,
  SAVE_ANGULAR_FLUX = 8,
  USE_SOURCE_MOMENTS = 9,
  VERBOSE_INNER_ITERATIONS = 10,
  VERBOSE_OUTER_ITERATIONS = 11,
  USE_PRECURSORS = 12,
};

int
chiLBSSetProperty(lua_State* L)
{
  const std::string fname = __FUNCTION__;

  opensn::Chi::log.Log0Warning() << fname + " has been deprecated. Use chiLBSSetOptions instead.";

  const int numArgs = lua_gettop(L);
  if (numArgs < 2) LuaPostArgAmountError(fname, 2, numArgs);

  LuaCheckNilValue(fname, L, 1);

  // Get pointer to solver
  const int solver_handle = lua_tonumber(L, 1);
  auto& lbs_solver = opensn::Chi::GetStackItem<opensn::lbs::LBSSolver>(
    opensn::Chi::object_stack, solver_handle, fname);

  // Get property index
  LuaCheckNilValue(fname, L, 2);

  const int property = lua_tonumber(L, 2);

  // Handle properties
  if (scpcode(property) == PropertyCode::DISCRETIZATION_METHOD)
  {
    LuaCheckNilValue(fname, L, 3);

    const int method = lua_tonumber(L, 3);

    typedef SpatialDiscretizationType SDMType;

    if (scpcode(method) == PropertyCode::PWLD)
      lbs_solver.Options().sd_type = SDMType::PIECEWISE_LINEAR_DISCONTINUOUS;
    else
      throw std::invalid_argument(
        "Invalid option for Discretization method in chiLBSSetProperty.\n");
  }
  else if (scpcode(property) == PropertyCode::BOUNDARY_CONDITION)
  {
    if (numArgs < 4) LuaPostArgAmountError("chiLBSSetProperty", 4, numArgs);

    LuaCheckNilValue(fname, L, 3);
    LuaCheckNilValue(fname, L, 4);

    const int bident = lua_tonumber(L, 3);
    const int btype = lua_tonumber(L, 4);

    if (!((bident >= scint(PropertyCode::XMAX)) && (bident <= scint(PropertyCode::ZMIN))))
    {
      opensn::Chi::log.LogAllError() << "Unknown boundary identifier encountered "
                                        "in call to chiLBSSetProperty";
      opensn::Chi::Exit(EXIT_FAILURE);
    }

    const int bid = bident - 31;

    if (btype == (int)opensn::lbs::BoundaryType::VACUUM)
    {
      lbs_solver.BoundaryPreferences()[bid] = {opensn::lbs::BoundaryType::VACUUM};
      opensn::Chi::log.Log() << "Boundary " << bid << " set to Vacuum.";
    }
    else if (btype == (int)opensn::lbs::BoundaryType::INCIDENT_ISOTROPIC)
    {
      if (numArgs != 5) LuaPostArgAmountError("chiLBSSetProperty", 5, numArgs);

      if (lbs_solver.Groups().empty())
      {
        opensn::Chi::log.Log0Error()
          << "In call to chiLBSSetProperty, setting "
          << "incident isotropic flux boundary type: Number of solver groups"
          << " is zero. Boundary fluxes can only be set after group structure"
          << " has been defined.";
        opensn::Chi::Exit(EXIT_FAILURE);
      }

      if (!lua_istable(L, 5))
      {
        opensn::Chi::log.LogAllError() << "In call to chiLBSSetProperty, setting "
                                       << "incident isotropic flux boundary type,"
                                       << " argument 5 should be a lua table and was detected as"
                                          " not being one.";
        opensn::Chi::Exit(EXIT_FAILURE);
      }

      const size_t table_len = lua_rawlen(L, 5);
      std::vector<double> group_strength(table_len, 0.0);
      for (int g = 0; g < table_len; g++)
      {
        lua_pushnumber(L, g + 1);
        lua_gettable(L, 5);
        group_strength[g] = lua_tonumber(L, -1);
        lua_pop(L, 1);
      }

      if (table_len != lbs_solver.Groups().size())
      {
        opensn::Chi::log.Log0Error()
          << "In call to chiLBSSetProperty, setting "
          << "incident isotropic flux boundary type: "
          << "Number of groups in boundary flux specification is " << table_len
          << " but solver has a total of " << lbs_solver.Groups().size()
          << " groups. These two must be equal.";
        opensn::Chi::Exit(EXIT_FAILURE);
      }

      lbs_solver.BoundaryPreferences()[bid] = {opensn::lbs::BoundaryType::INCIDENT_ISOTROPIC,
                                               group_strength};

      opensn::Chi::log.Log() << "Isotropic boundary condition for boundary " << bid
                             << " loaded with " << table_len << " groups.";
    }
    else if (btype == (int)opensn::lbs::BoundaryType::REFLECTING)
    {
      lbs_solver.BoundaryPreferences()[bid] = {opensn::lbs::BoundaryType::REFLECTING};
      opensn::Chi::log.Log() << "Boundary " << bid << " set to Reflecting.";
    }
    else if (btype == (int)opensn::lbs::BoundaryType::INCIDENT_ANISTROPIC_HETEROGENEOUS)
    {
      LuaCheckNilValue(fname, L, 5);

      const std::string lua_func_name = lua_tostring(L, 5);
      lbs_solver.BoundaryPreferences()[bid] = {
        opensn::lbs::BoundaryType::INCIDENT_ANISTROPIC_HETEROGENEOUS, {}, lua_func_name};
      opensn::Chi::log.Log() << "Boundary " << bid
                             << " set to Incident anistoropic"
                                " heterogeneous.";
    }
    else
    {
      opensn::Chi::log.LogAllError() << "Unsupported boundary type encountered "
                                        "in call to "
                                     << LuaSourceInfo(L, "chiLBSSetProperty");
      opensn::Chi::Exit(EXIT_FAILURE);
    }
  }
  else if (scpcode(property) == PropertyCode::SCATTERING_ORDER)
  {
    LuaCheckNilValue(fname, L, 3);

    const int scattering_order = lua_tonumber(L, 3);

    if (scattering_order < 0)
    {
      opensn::Chi::log.Log0Error() << "Invalid scattering order in call to "
                                   << "chiLBSSetProperty:SCATTERING_ORDER. "
                                      "Value must be > 0.";
      opensn::Chi::Exit(EXIT_FAILURE);
    }

    lbs_solver.Options().scattering_order = scattering_order;
  }
  else if (scpcode(property) == PropertyCode::SWEEP_EAGER_LIMIT)
  {
    if (numArgs != 3) LuaPostArgAmountError("chiLBSSetProperty:SWEEP_EAGER_LIMIT", 3, numArgs);

    LuaCheckNilValue(fname, L, 3);

    const int limit = lua_tonumber(L, 3);
    lbs_solver.Options().sweep_eager_limit = limit;
  }
  else if (scpcode(property) == PropertyCode::READ_RESTART_DATA)
  {
    if (numArgs >= 3)
    {
      LuaCheckNilValue(fname, L, 3);

      const std::string folder = lua_tostring(L, 3);
      lbs_solver.Options().read_restart_folder_name = std::string(folder);
      opensn::Chi::log.Log() << "Restart input folder set to " << folder;
    }
    if (numArgs >= 4)
    {
      LuaCheckNilValue(fname, L, 4);

      const std::string filebase = lua_tostring(L, 4);
      lbs_solver.Options().read_restart_file_base = std::string(filebase);
      opensn::Chi::log.Log() << "Restart input filebase set to " << filebase;
    }
    lbs_solver.Options().read_restart_data = true;
  }
  else if (scpcode(property) == PropertyCode::WRITE_RESTART_DATA)
  {
    if (numArgs >= 3)
    {
      LuaCheckNilValue(fname, L, 3);

      const std::string folder = lua_tostring(L, 3);
      lbs_solver.Options().write_restart_folder_name = std::string(folder);
      opensn::Chi::log.Log() << "Restart output folder set to " << folder;
    }
    if (numArgs >= 4)
    {
      LuaCheckNilValue(fname, L, 4);

      const std::string filebase = lua_tostring(L, 4);
      lbs_solver.Options().write_restart_file_base = std::string(filebase);
      opensn::Chi::log.Log() << "Restart output filebase set to " << filebase;
    }
    if (numArgs == 5)
    {
      LuaCheckNilValue(fname, L, 5);

      const double interval = lua_tonumber(L, 5);
      lbs_solver.Options().write_restart_interval = interval;
    }
    lbs_solver.Options().write_restart_data = true;
  }
  else if (scpcode(property) == PropertyCode::SAVE_ANGULAR_FLUX)
  {
    LuaCheckNilValue(fname, L, 3);

    const bool save_flag = lua_toboolean(L, 3);

    lbs_solver.Options().save_angular_flux = save_flag;

    opensn::Chi::log.Log() << "LBS option to save angular flux set to " << save_flag;
  }
  else if (scpcode(property) == PropertyCode::USE_SOURCE_MOMENTS)
  {
    LuaCheckNilValue(fname, L, 3);

    const bool use_flag = lua_toboolean(L, 3);

    lbs_solver.Options().use_src_moments = use_flag;

    opensn::Chi::log.Log() << "LBS option to use source moments set to " << use_flag;
  }
  else if (scpcode(property) == PropertyCode::VERBOSE_INNER_ITERATIONS)
  {
    LuaCheckNilValue(fname, L, 3);

    const bool flag = lua_toboolean(L, 3);

    lbs_solver.Options().verbose_inner_iterations = flag;

    opensn::Chi::log.Log() << "LBS option: verbose_inner_iterations set to " << flag;
  }
  else if (scpcode(property) == PropertyCode::VERBOSE_OUTER_ITERATIONS)
  {
    LuaCheckNilValue(fname, L, 3);

    const bool flag = lua_toboolean(L, 3);

    lbs_solver.Options().verbose_outer_iterations = flag;

    opensn::Chi::log.Log() << "LBS option: verbose_outer_iterations set to " << flag;
  }
  else if (scpcode(property) == PropertyCode::USE_PRECURSORS)
  {
    LuaCheckNilValue(fname, L, 3);

    const bool flag = lua_toboolean(L, 3);

    lbs_solver.Options().use_precursors = flag;

    opensn::Chi::log.Log() << "LBS option: use_precursors set to " << flag;
  }
  else
    throw std::logic_error(fname + ": Invalid property in chiLBSSetProperty.\n");

  return 0;
}

} // namespace opensnlua::lbs
