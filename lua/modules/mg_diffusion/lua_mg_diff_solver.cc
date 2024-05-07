// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "lua_mg_diff_solver.h"
#include "lua/framework/console/console.h"
#include "modules/mg_diffusion/mg_diffusion_solver.h"
#include "framework/runtime.h"
#include "framework/logging/log.h"

using namespace opensn;

namespace opensnlua
{

RegisterLuaFunctionInNamespace(CFEMMGDiffusionSolverCreate, diffusion, CFEMMGSolverCreate);
RegisterLuaFunctionInNamespace(CFEMMGDiffusionSetBCProperty, diffusion, CFEMMGSetBCProperty);

int
CFEMMGDiffusionSolverCreate(lua_State* L)
{
  const std::string fname = "diffusion.CFEMMGSolverCreate";

  auto solver_name = LuaArgOptional<std::string>(L, 1, "MGDiffusionSolver");

  auto new_solver = std::make_shared<opensn::mg_diffusion::Solver>(solver_name);

  opensn::object_stack.push_back(new_solver);

  opensn::log.LogAllVerbose1() << fname << ": CFEM Multigroup Diffusion solver created.";

  return LuaReturn(L, opensn::object_stack.size() - 1);
}

int
CFEMMGDiffusionSetBCProperty(lua_State* L)
{
  const std::string fname = "diffusion.CFEMMGSetBCProperty";
  LuaCheckArgs<size_t, std::string>(L, fname);

  // Get solver
  const auto solver_index = LuaArg<size_t>(L, 1);
  auto& solver =
    opensn::GetStackItem<opensn::mg_diffusion::Solver>(opensn::object_stack, solver_index, fname);

  // Get property index
  const auto property_name = LuaArg<std::string>(L, 2);

  // Handle properties
  if (property_name == "boundary_type")
  {
    const auto bound_index = LuaArg<int>(L, 3);
    const auto type_name = LuaArg<std::string>(L, 4);

    if (type_name == "reflecting") // ------------- REFLECTING
    {
      opensn::mg_diffusion::Solver::BoundaryInfo bndry_info;
      bndry_info.first = opensn::mg_diffusion::BoundaryType::Reflecting;

      solver.boundary_preferences_.insert(std::make_pair(bound_index, bndry_info));

      opensn::log.Log() << "Boundary " << bound_index << " set as "
                        << "Reflecting.";
    }
    else if (type_name == "dirichlet") // ------------- DIRICHLET
    {
      opensn::log.Log0Error() << fname
                              << ": Dirichlet BC is not supported in multigroup diffusion.";
      opensn::Exit(EXIT_FAILURE);
    }
    else if (type_name == "neumann") // ------------- NEUMANN
    {
      auto f_values = LuaArg<std::vector<double>>(L, 5);
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
      // check lua tables
      auto a_values = LuaArg<std::vector<double>>(L, 5);
      auto b_values = LuaArg<std::vector<double>>(L, 6);
      auto f_values = LuaArg<std::vector<double>>(L, 7);

      opensn::mg_diffusion::Solver::BoundaryInfo bndry_info;
      bndry_info.first = opensn::mg_diffusion::BoundaryType::Robin;
      bndry_info.second = {a_values, b_values, f_values};
      solver.boundary_preferences_.insert(std::make_pair(bound_index, bndry_info));

      opensn::log.Log() << "Boundary " << bound_index << " set as Robin";
    }
    else
    {
      opensn::log.LogAllError() << fname << ": Unsupported boundary type '" << type_name << "'.";
      opensn::Exit(EXIT_FAILURE);
    }
  }
  else // wrong property_name
  {
    opensn::log.Log0Error() << fname + ": Invalid property '" + property_name + "'.";
    opensn::Exit(EXIT_FAILURE);
  }
  return LuaReturn(L);
}

} // namespace opensnlua
