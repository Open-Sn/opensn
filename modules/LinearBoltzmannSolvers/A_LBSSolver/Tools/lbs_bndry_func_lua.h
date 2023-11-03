#pragma once

#include "mesh/SweepUtilities/SweepBoundary/sweep_boundaries.h"

#include <string>
#include <utility>

namespace lbs
{

class BoundaryFunctionToLua : public chi_mesh::sweep_management::BoundaryFunction
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
           const chi_mesh::Vector3& face_node_location,
           const chi_mesh::Vector3& face_node_normal,
           const std::vector<int>& quadrature_angle_indices,
           const std::vector<chi_mesh::Vector3>& quadrature_angle_vectors,
           const std::vector<std::pair<double, double>>& quadrature_phi_theta_angles,
           const std::vector<int>& group_indices,
           double time) override;
};

} // namespace lbs


