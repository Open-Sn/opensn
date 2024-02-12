#include "framework/lua.h"

#include "modules/mg_diffusion/mg_diffusion_solver.h"

#include "framework/runtime.h"
#include "framework/logging/log.h"

using namespace opensn;

namespace opensnlua::mg_diffusion
{

int
CFEMMGDiffusionSetBCProperty(lua_State* L)
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
    opensn::GetStackItem<opensn::mg_diffusion::Solver>(opensn::object_stack, solver_index, fname);

  // Get property index
  LuaCheckStringValue(fname, L, 2);
  const std::string property_name = lua_tostring(L, 2);

  // Handle properties
  if (property_name == "boundary_type")
  {
    if (num_args < 4)
    {
      opensn::log.Log0Error() << "Invalid amount of arguments used in"
                              << " CFEMMGDiffusionsetBCproperty(...,\"boundary_type\".... "
                              << " At least 4 arguments are expected.";
      opensn::Exit(EXIT_FAILURE);
    }
    LuaCheckNumberValue(fname, L, 3);
    const int bound_index = lua_tonumber(L, 3);

    LuaCheckStringValue(fname, L, 4);
    const std::string type_name = lua_tostring(L, 4);

    if (type_name == "reflecting") // ------------- REFLECTING
    {
      if (num_args != 4)
      {
        opensn::log.Log0Error() << "Invalid amount of arguments used in"
                                << " CFEMMGDiffusionsetBCproperty(...,\"boundary_type\","
                                << bound_index << ",\"reflecting\". "
                                << " 4 arguments are expected.";
        opensn::Exit(EXIT_FAILURE);
      }

      opensn::mg_diffusion::Solver::BoundaryInfo bndry_info;
      bndry_info.first = opensn::mg_diffusion::BoundaryType::Reflecting;

      solver.boundary_preferences_.insert(std::make_pair(bound_index, bndry_info));

      opensn::log.Log() << "Boundary " << bound_index << " set as "
                        << "Reflecting.";
    }
    else if (type_name == "dirichlet") // ------------- DIRICHLET
    {
      opensn::log.Log0Error() << "Dirichlet BC is not supported in multigroup diffusion "
                              << "(CFEMMGDiffusionSetBCProperty).";
      opensn::Exit(EXIT_FAILURE);
    }
    else if (type_name == "neumann") // ------------- NEUMANN
    {
      if (num_args != 5)
      {
        opensn::log.Log0Error() << "Invalid amount of arguments used in"
                                << " CFEMMGDiffusionsetBCproperty(...,\"boundary_type\","
                                << bound_index << ",\"neumann\". "
                                << " 5 arguments are expected.";
        opensn::Exit(EXIT_FAILURE);
      }
      // check lua tables
      LuaCheckTableValue(fname, L, 5);
      std::vector<double> f_values;
      LuaPopulateVectorFrom1DArray(fname, L, 5, f_values);
      // add the other multigroup vectors to finish the BC
      unsigned int ng = f_values.size();
      std::vector<double> a_values(ng, 0.0);
      std::vector<double> b_values(ng, 1.0);

      opensn::mg_diffusion::Solver::BoundaryInfo bndry_info;
      bndry_info.first = opensn::mg_diffusion::BoundaryType::Neumann;
      bndry_info.second = {a_values, b_values, f_values};
      solver.boundary_preferences_.insert(std::make_pair(bound_index, bndry_info));

      opensn::log.Log() << "Boundary " << bound_index << " set as "
                        << "Neumann with D_g grad(u_g) dot n = f_g";
    }
    else if (type_name == "vacuum") // ------------- VACUUM
    {
      if (num_args != 4)
      {
        opensn::log.Log0Error() << "Invalid amount of arguments used in"
                                << " CFEMMGDiffusionsetBCproperty(...,\"boundary_type\","
                                << bound_index << ",\"vacuum\". "
                                << " 4 arguments are expected.";
        opensn::Exit(EXIT_FAILURE);
      }
      // dummy-sized values until we now num_group later, after solver init
      std::vector<double> a_values(1, 0.25);
      std::vector<double> b_values(1, 0.5);
      std::vector<double> f_values(1, 0.0);

      opensn::mg_diffusion::Solver::BoundaryInfo bndry_info;
      bndry_info.first = opensn::mg_diffusion::BoundaryType::Vacuum;
      bndry_info.second = {a_values, b_values, f_values};
      solver.boundary_preferences_.insert(std::make_pair(bound_index, bndry_info));

      opensn::log.Log() << "Boundary " << bound_index << " set as "
                        << "Vacuum.";
    }
    else if (type_name == "robin") // ------------- ROBIN
    {
      if (num_args != 7)
      {
        opensn::log.Log0Error() << "Invalid amount of arguments used in"
                                << " CFEMMGDiffusionSetBCProperty(...,\"boundary_type\","
                                << bound_index << ",\"robin\". "
                                << " 7 arguments are expected.";
        opensn::Exit(EXIT_FAILURE);
      }
      // check lua tables
      LuaCheckTableValue(fname, L, 5);
      LuaCheckTableValue(fname, L, 6);
      LuaCheckTableValue(fname, L, 7);
      std::vector<double> a_values, b_values, f_values;
      LuaPopulateVectorFrom1DArray(fname, L, 5, a_values);
      LuaPopulateVectorFrom1DArray(fname, L, 6, b_values);
      LuaPopulateVectorFrom1DArray(fname, L, 7, f_values);

      opensn::mg_diffusion::Solver::BoundaryInfo bndry_info;
      bndry_info.first = opensn::mg_diffusion::BoundaryType::Robin;
      bndry_info.second = {a_values, b_values, f_values};
      solver.boundary_preferences_.insert(std::make_pair(bound_index, bndry_info));

      opensn::log.Log() << "Boundary " << bound_index << " set as Robin";
    }
    else
    {
      opensn::log.LogAllError() << "Unsupported boundary type encountered in call to "
                                << "CFEMMGDiffusionSetBCProperty(..,\"boundary_type\",.. :"
                                << type_name;
      opensn::Exit(EXIT_FAILURE);
    }
  }
  else // wrong property_name
  {
    opensn::log.Log0Error() << "Invalid property in CFEMMGDiffusionSetBCProperty.";
    opensn::Exit(EXIT_FAILURE);
  }
  return 0;
} // end of CFEMMGDiffusionSetBCProperty

} // namespace opensnlua::mg_diffusion
