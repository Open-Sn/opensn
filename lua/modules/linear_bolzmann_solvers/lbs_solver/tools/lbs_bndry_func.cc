// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "lua/modules/linear_bolzmann_solvers/lbs_solver/tools/lbs_bndry_func.h"
#include "lua/framework/lua.h"
#include "framework/runtime.h"
#include "framework/logging/log.h"
#include "lua/framework/console/console.h"

using namespace opensn;

namespace opensnlua
{
namespace lbs
{

std::vector<double>
BoundaryFunctionToLua::Evaluate(
  size_t cell_global_id,
  int cell_material_id,
  unsigned int face_index,
  unsigned int face_node_index,
  const Vector3& face_node_location,
  const Vector3& face_node_normal,
  const std::vector<int>& quadrature_angle_indices,
  const std::vector<Vector3>& quadrature_angle_vectors,
  const std::vector<std::pair<double, double>>& quadrature_phi_theta_angles,
  const std::vector<int>& group_indices,
  double time)
{
  const std::string fname = "LinearBoltzmann::BoundaryFunctionToLua";

  // Get lua function
  lua_State* L = console.GetConsoleState();
  auto psi = LuaCall<std::vector<double>>(L,
                                          m_lua_function_name,
                                          cell_global_id,
                                          cell_material_id,
                                          face_node_location,
                                          face_node_normal,
                                          quadrature_angle_indices,
                                          quadrature_angle_vectors,
                                          quadrature_phi_theta_angles,
                                          group_indices,
                                          time);

  // Error check psi vector
  size_t num_angles = quadrature_angle_indices.size();
  size_t num_groups = group_indices.size();

  if (psi.size() != (num_angles * num_groups))
    throw std::logic_error(fname + " the returned vector from lua-function, " +
                           m_lua_function_name + ", did not produce the required size vector. " +
                           "The size must equal num_angles*num_groups, " +
                           std::to_string(num_angles * num_groups) + ", but the size is " +
                           std::to_string(psi.size()) + ".");

  return psi;
}

} // namespace lbs
} // namespace opensnlua
