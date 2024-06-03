// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "lua/modules/fv_diffusion/fv_diff_solver.h"
#include "modules/fv_diffusion/fv_diffusion_solver.h"
#include "lua/framework/console/console.h"
#include "framework/runtime.h"
#include "framework/logging/log.h"
#include "lua/framework/math/functions/lua_scalar_spatial_material_function.h"

using namespace opensn;

namespace opensnlua
{

RegisterLuaFunctionInNamespace(FVDiffusionSolverCreate, diffusion, FVSolverCreate);
RegisterLuaFunctionInNamespace(FVDiffusionSetBCProperty, diffusion, FVSetBCProperty);

namespace
{

std::shared_ptr<LuaScalarSpatialMaterialFunction>
CreateFunction(const std::string& function_name)
{
  ParameterBlock blk;
  blk.AddParameter("lua_function_name", function_name);
  InputParameters params = LuaScalarSpatialMaterialFunction::GetInputParameters();
  params.AssignParameters(blk);
  return std::make_shared<LuaScalarSpatialMaterialFunction>(params);
}

} // namespace

int
FVDiffusionSolverCreate(lua_State* L)
{
  const std::string fname = __FUNCTION__;

  auto solver_name = LuaArgOptional<std::string>(L, 1, "FVDiffusionSolver");

  auto d_coef_function = CreateFunction("D_coef");
  opensn::function_stack.push_back(d_coef_function);

  auto q_ext_function = CreateFunction("Q_ext");
  opensn::function_stack.push_back(q_ext_function);

  auto sigma_a_function = CreateFunction("Sigma_a");
  opensn::function_stack.push_back(sigma_a_function);

  auto new_solver = std::make_shared<opensn::fv_diffusion::Solver>(solver_name);
  new_solver->SetDCoefFunction(d_coef_function);
  new_solver->SetQExtFunction(q_ext_function);
  new_solver->SetSigmaAFunction(sigma_a_function);

  opensn::object_stack.push_back(new_solver);

  opensn::log.LogAllVerbose1() << "\nFVDiffusionSolverCreate: FV Diffusion solver created"
                               << std::endl;

  return LuaReturn(L, opensn::object_stack.size() - 1);
}

int
FVDiffusionSetBCProperty(lua_State* L)
{
  const std::string fname = "diffusion.FVSetBCProperty";
  LuaCheckArgs<size_t>(L, fname);

  // Get solver
  const auto solver_index = LuaArg<size_t>(L, 1);
  auto& solver =
    opensn::GetStackItem<opensn::fv_diffusion::Solver>(opensn::object_stack, solver_index, fname);

  // Get property index
  const auto property_name = LuaArg<std::string>(L, 2);

  // Handle properties
  if (property_name == "boundary_type")
  {
    const auto bound_name = LuaArg<std::string>(L, 3);
    const auto type_name = LuaArg<std::string>(L, 4);

    if (type_name == "reflecting")
    {
      opensn::fv_diffusion::Solver::BoundaryInfo bndry_info;
      bndry_info.first = opensn::fv_diffusion::BoundaryType::Reflecting;

      solver.boundary_preferences_.insert(std::make_pair(bound_name, bndry_info));

      opensn::log.Log() << "Boundary " << bound_name << " set as "
                        << "Reflecting.";
    }
    else if (type_name == "dirichlet")
    {
      auto boundary_value = LuaArg<double>(L, 5);

      opensn::fv_diffusion::Solver::BoundaryInfo bndry_info;
      bndry_info.first = opensn::fv_diffusion::BoundaryType::Dirichlet;
      bndry_info.second = {boundary_value};
      solver.boundary_preferences_.insert(std::make_pair(bound_name, bndry_info));

      opensn::log.Log() << "Boundary " << bound_name << " set as "
                        << "Dirichlet with value " << boundary_value;
    }
    else if (type_name == "neumann")
    {
      auto f_value = LuaArg<double>(L, 5);

      opensn::fv_diffusion::Solver::BoundaryInfo bndry_info;
      bndry_info.first = opensn::fv_diffusion::BoundaryType::Robin;
      bndry_info.second = {0.0, 1.0, f_value};
      solver.boundary_preferences_.insert(std::make_pair(bound_name, bndry_info));

      opensn::log.Log() << "Boundary " << bound_name << " set as "
                        << "Neumann with D grad(u) dot n = (" << f_value << ") ";
    }
    else if (type_name == "vacuum")
    {
      opensn::fv_diffusion::Solver::BoundaryInfo bndry_info;
      bndry_info.first = opensn::fv_diffusion::BoundaryType::Robin;
      bndry_info.second = {0.25, 0.5, 0.0};
      solver.boundary_preferences_.insert(std::make_pair(bound_name, bndry_info));

      opensn::log.Log() << "Boundary " << bound_name << " set as "
                        << "Vacuum.";
    }
    else if (type_name == "robin")
    {
      auto a_value = LuaArg<double>(L, 5);
      auto b_value = LuaArg<double>(L, 6);
      auto f_value = LuaArg<double>(L, 7);

      opensn::fv_diffusion::Solver::BoundaryInfo bndry_info;
      bndry_info.first = opensn::fv_diffusion::BoundaryType::Robin;
      bndry_info.second = {a_value, b_value, f_value};
      solver.boundary_preferences_.insert(std::make_pair(bound_name, bndry_info));

      opensn::log.Log() << "Boundary " << bound_name << " set as "
                        << "Robin with a,b,f = (" << a_value << "," << b_value << "," << f_value
                        << ") ";
    }
    else
    {
      opensn::log.LogAllError() << fname << ": Unsupported boundary type '" << type_name << "'.";
      opensn::Exit(EXIT_FAILURE);
    }
  }
  else
  {
    opensn::log.Log0Error() << fname + ": Invalid property '" + property_name + "'.";
    opensn::Exit(EXIT_FAILURE);
  }

  return LuaReturn(L);
}

} // namespace opensnlua
