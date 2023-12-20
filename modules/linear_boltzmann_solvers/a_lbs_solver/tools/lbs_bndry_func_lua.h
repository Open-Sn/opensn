#pragma once
#ifdef OPENSN_WITH_LUA
#include "framework/mesh/sweep_utilities/sweep_boundary/sweep_boundary.h"

#include <string>
#include <utility>

namespace opensn
{
namespace lbs
{

class BoundaryFunctionToLua : public BoundaryFunction
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
           const Vector3& face_node_location,
           const Vector3& face_node_normal,
           const std::vector<int>& quadrature_angle_indices,
           const std::vector<Vector3>& quadrature_angle_vectors,
           const std::vector<std::pair<double, double>>& quadrature_phi_theta_angles,
           const std::vector<int>& group_indices,
           double time) override;
};

} // namespace lbs
} // namespace opensn
#endif
