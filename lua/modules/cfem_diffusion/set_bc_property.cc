#include "framework/lua.h"

#include "modules/cfem_diffusion/cfem_diffusion_solver.h"

#include "framework/runtime.h"
#include "framework/logging/log.h"

using namespace opensn;

namespace opensnlua::cfem_diffusion
{

int
CFEMDiffusionSetBCProperty(lua_State* L)
{
  const std::string fname = __FUNCTION__;
  const int num_args = lua_gettop(L);
  if (num_args < 2)
    LuaPostArgAmountError(fname, num_args, 2);

  LuaCheckNilValue(fname, L, 1);
  LuaCheckNilValue(fname, L, 2);

  // Get solver
  LuaCheckNumberValue(fname, L, 1);
  const int solver_index = lua_tonumber(L, 1);

  auto& solver =
    opensn::GetStackItem<opensn::cfem_diffusion::Solver>(opensn::object_stack, solver_index, fname);

  // Get property index
  LuaCheckStringValue(fname, L, 2);
  const std::string property_name = lua_tostring(L, 2);

  // Handle properties
  if (property_name == "boundary_type")
  {
    if (num_args < 4)
    {
      opensn::log.Log0Error() << "Invalid amount of arguments used in"
                              << " CFEMDiffusionsetBCproperty(...,\"boundary_type\".... "
                              << " At least 4 arguments are expected.";
      opensn::Exit(EXIT_FAILURE);
    }
    LuaCheckStringValue(fname, L, 3);
    const std::string bound_name = lua_tostring(L, 3);

    LuaCheckStringValue(fname, L, 4);
    const std::string type_name = lua_tostring(L, 4);

    if (type_name == "reflecting")
    {
      if (num_args != 4)
      {
        opensn::log.Log0Error() << "Invalid amount of arguments used in"
                                << " CFEMDiffusionsetBCproperty(...,\"boundary_type\","
                                << bound_name << ",\"reflecting\". "
                                << " 4 arguments are expected.";
        opensn::Exit(EXIT_FAILURE);
      }

      opensn::cfem_diffusion::Solver::BoundaryInfo bndry_info;
      bndry_info.first = opensn::cfem_diffusion::BoundaryType::Reflecting;

      solver.boundary_preferences_.insert(std::make_pair(bound_name, bndry_info));

      opensn::log.Log() << "Boundary " << bound_name << " set as "
                        << "Reflecting.";
    }
    else if (type_name == "dirichlet")
    {
      if (num_args != 5)
      {
        opensn::log.Log0Error() << "Invalid amount of arguments used in"
                                << " CFEMDiffusionsetBCproperty(...,\"boundary_type\","
                                << bound_name << ",\"dirichlet\". "
                                << " 5 arguments are expected.";
        opensn::Exit(EXIT_FAILURE);
      }
      LuaCheckNumberValue(fname, L, 5);
      double boundary_value = lua_tonumber(L, 5);

      opensn::cfem_diffusion::Solver::BoundaryInfo bndry_info;
      bndry_info.first = opensn::cfem_diffusion::BoundaryType::Dirichlet;
      bndry_info.second = {boundary_value};
      solver.boundary_preferences_.insert(std::make_pair(bound_name, bndry_info));

      opensn::log.Log() << "Boundary " << bound_name << " set as "
                        << "Dirichlet with value " << boundary_value;
    }
    else if (type_name == "neumann")
    {
      if (num_args != 5)
      {
        opensn::log.Log0Error() << "Invalid amount of arguments used in"
                                << " CFEMDiffusionsetBCproperty(...,\"boundary_type\","
                                << bound_name << ",\"neumann\". "
                                << " 5 arguments are expected.";
        opensn::Exit(EXIT_FAILURE);
      }
      LuaCheckNumberValue(fname, L, 5);
      double f_value = lua_tonumber(L, 5);

      opensn::cfem_diffusion::Solver::BoundaryInfo bndry_info;
      bndry_info.first = opensn::cfem_diffusion::BoundaryType::Robin;
      bndry_info.second = {0.0, 1.0, f_value};
      solver.boundary_preferences_.insert(std::make_pair(bound_name, bndry_info));

      opensn::log.Log() << "Boundary " << bound_name << " set as "
                        << "Neumann with D grad(u) dot n = (" << f_value << ") ";
    }
    else if (type_name == "vacuum")
    {
      if (num_args != 4)
      {
        opensn::log.Log0Error() << "Invalid amount of arguments used in"
                                << " CFEMDiffusionsetBCproperty(...,\"boundary_type\","
                                << bound_name << ",\"vacuum\". "
                                << " 4 arguments are expected.";
        opensn::Exit(EXIT_FAILURE);
      }

      opensn::cfem_diffusion::Solver::BoundaryInfo bndry_info;
      bndry_info.first = opensn::cfem_diffusion::BoundaryType::Robin;
      bndry_info.second = {0.25, 0.5, 0.0};
      solver.boundary_preferences_.insert(std::make_pair(bound_name, bndry_info));

      opensn::log.Log() << "Boundary " << bound_name << " set as "
                        << "Vacuum.";
    }
    else if (type_name == "robin")
    {
      if (num_args != 7)
      {
        opensn::log.Log0Error() << "Invalid amount of arguments used in"
                                << " CFEMDiffusionsetBCproperty(...,\"boundary_type\","
                                << bound_name << ",\"robin\". "
                                << " 7 arguments are expected.";
        opensn::Exit(EXIT_FAILURE);
      }
      LuaCheckNumberValue(fname, L, 5);
      LuaCheckNumberValue(fname, L, 6);
      LuaCheckNumberValue(fname, L, 7);

      double a_value = lua_tonumber(L, 5);
      double b_value = lua_tonumber(L, 6);
      double f_value = lua_tonumber(L, 7);

      opensn::cfem_diffusion::Solver::BoundaryInfo bndry_info;
      bndry_info.first = opensn::cfem_diffusion::BoundaryType::Robin;
      bndry_info.second = {a_value, b_value, f_value};
      solver.boundary_preferences_.insert(std::make_pair(bound_name, bndry_info));

      opensn::log.Log() << "Boundary " << bound_name << " set as "
                        << "Robin with a,b,f = (" << a_value << "," << b_value << "," << f_value
                        << ") ";
    }
    else
    {
      opensn::log.LogAllError() << "Unsupported boundary type encountered in call to "
                                << "CFEMDiffusionSetBCProperty(..,\"boundary_type\",.. :"
                                << type_name;
      opensn::Exit(EXIT_FAILURE);
    }
  }
  else
  {
    opensn::log.Log0Error() << "Invalid property in DiffusionsetBCproperty.";
    opensn::Exit(EXIT_FAILURE);
  }
  return 0;
}

} // namespace opensnlua::cfem_diffusion
