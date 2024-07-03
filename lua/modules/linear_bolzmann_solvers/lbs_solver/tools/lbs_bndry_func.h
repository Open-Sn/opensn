// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include "modules/linear_boltzmann_solvers/discrete_ordinates_solver/sweep/boundary/sweep_boundary.h"
#include <string>
#include <utility>

namespace opensnlua
{

class BoundaryFunctionToLua : public opensn::BoundaryFunction
{
private:
  const std::string m_lua_function_name;

public:
  explicit BoundaryFunctionToLua(std::string lua_function_name)
    : m_lua_function_name(std::move(lua_function_name))
  {
  }

  std::vector<double>
  Evaluate(size_t cell_global_id,
           int cell_material_id,
           unsigned int face_index,
           unsigned int face_node_index,
           const opensn::Vector3& face_node_location,
           const opensn::Vector3& face_node_normal,
           const std::vector<int>& quadrature_angle_indices,
           const std::vector<opensn::Vector3>& quadrature_angle_vectors,
           const std::vector<std::pair<double, double>>& quadrature_phi_theta_angles,
           const std::vector<int>& group_indices,
           double time) override;
};

} // namespace opensnlua
