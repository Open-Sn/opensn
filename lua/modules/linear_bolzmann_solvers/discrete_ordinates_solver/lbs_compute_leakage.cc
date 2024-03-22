#include "lbs_do_lua_utils.h"
#include "framework/runtime.h"
#include "framework/console/console.h"
#include "framework/mesh/mesh_continuum/mesh_continuum.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_solver/lbs_discrete_ordinates_solver.h"

namespace opensnlua::lbs
{

RegisterLuaFunctionNamespace(ComputeLeakage, lbs, ComputeLeakage);

int
ComputeLeakage(lua_State* L)
{
  const auto fname = "lbs.ComputeLeakage";
  const auto num_args = lua_gettop(L);

  // Get the solver
  const auto solver_handle = LuaArg<size_t>(L, 1);
  const auto& solver = opensn::GetStackItem<opensn::lbs::DiscreteOrdinatesSolver>(
    opensn::object_stack, solver_handle, fname);

  // Get the supported boundaries
  const auto supported_boundary_names =
    opensn::lbs::DiscreteOrdinatesSolver::supported_boundary_names;
  const auto supported_boundary_ids = opensn::lbs::DiscreteOrdinatesSolver::supported_boundary_ids;

  // Get the boundaries to parse
  std::vector<uint64_t> bndry_ids;
  if (num_args > 1)
  {
    auto bnd_names = LuaArg<std::vector<std::string>>(L, 2);
    for (auto& name : bnd_names)
      bndry_ids.push_back(supported_boundary_names.at(name));
  }
  else
    bndry_ids = solver.Grid().GetDomainUniqueBoundaryIDs();

  // Compute the leakage
  const auto leakage = solver.ComputeLeakage(bndry_ids);

  // Push to lua table
  lua_newtable(L);
  for (const auto& [bid, vals] : leakage)
  {
    LuaPush(L, supported_boundary_ids.at(bid));

    lua_newtable(L);
    for (int g = 0; g < solver.NumGroups(); ++g)
      LuaPushTableKey(L, g + 1, vals[g]);
    lua_settable(L, -3);
  }
  return 1;
}

} // namespace opensnlua::lbs
