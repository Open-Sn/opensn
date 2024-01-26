#include "framework/lua.h"

#include "modules/cfem_diffusion/cfem_diffusion_solver.h"
#include "framework/runtime.h"
#include "framework/logging/log.h"
#include "lua/framework/math/functions/lua_scalar_spatial_material_function.h"

using namespace opensn;

namespace opensnlua::cfem_diffusion
{

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
chiCFEMDiffusionSolverCreate(lua_State* L)
{
  const std::string fname = __FUNCTION__;
  int num_args = lua_gettop(L);

  std::string solver_name = "CFEMDiffusionSolver";

  if (num_args == 1)
  {
    LuaCheckStringValue(fname, L, 1);
    solver_name = lua_tostring(L, 1);
  }

  auto d_coef_function = CreateFunction("D_coef");
  opensn::function_stack.push_back(d_coef_function);

  auto q_ext_function = CreateFunction("Q_ext");
  opensn::function_stack.push_back(q_ext_function);

  auto sigma_a_function = CreateFunction("Sigma_a");
  opensn::function_stack.push_back(sigma_a_function);

  auto new_solver = std::make_shared<opensn::cfem_diffusion::Solver>(solver_name);
  new_solver->SetDCoefFunction(d_coef_function);
  new_solver->SetQExtFunction(q_ext_function);
  new_solver->SetSigmaAFunction(sigma_a_function);

  opensn::object_stack.push_back(new_solver);

  lua_pushinteger(L, static_cast<lua_Integer>(opensn::object_stack.size() - 1));

  opensn::log.LogAllVerbose1() << "\nCFEMDiffusionSolverCreate: CFEM Diffusion solver created"
                               << std::endl;
  return 1;
}

} // namespace opensnlua::cfem_diffusion
