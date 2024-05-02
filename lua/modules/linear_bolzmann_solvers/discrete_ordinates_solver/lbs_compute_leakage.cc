// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "lbs_do_lua_utils.h"
#include "framework/runtime.h"
#include "framework/console/console.h"
#include "framework/mesh/mesh_continuum/mesh_continuum.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_solver/lbs_discrete_ordinates_solver.h"

namespace opensnlua::lbs
{

RegisterLuaFunctionInNamespace(ComputeLeakage, lbs, ComputeLeakage);

int
ComputeLeakage(lua_State* L)
{
  const auto fname = "lbs.ComputeLeakage";
  LuaCheckArgs<size_t>(L, fname);

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
  if (LuaNumArgs(L) > 1)
  {
    auto bnd_names = LuaArg<std::vector<std::string>>(L, 2);
    for (auto& name : bnd_names)
      bndry_ids.push_back(supported_boundary_names.at(name));
  }
  else
    bndry_ids = solver.Grid().GetDomainUniqueBoundaryIDs();

  // Compute the leakage
  const auto leakage = solver.ComputeLeakage(bndry_ids);

  std::map<std::string, std::vector<double>> ret_val;
  for (const auto& [bid, vals] : leakage)
  {
    auto bnd_name = supported_boundary_ids.at(bid);
    ret_val.insert(std::pair<std::string, std::vector<double>>(bnd_name, vals));
  }
  return LuaReturn(L, ret_val);
}

} // namespace opensnlua::lbs
