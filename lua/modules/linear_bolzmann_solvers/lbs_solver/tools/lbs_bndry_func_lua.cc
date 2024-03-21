#include "lua/modules/linear_bolzmann_solvers/lbs_solver/tools/lbs_bndry_func_lua.h"
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

  auto PushPhiThetaPairTable = [](lua_State* L, const std::pair<double, double>& phi_theta)
  {
    lua_newtable(L);
    LuaPushTableKey(L, "phi", phi_theta.first);
    LuaPushTableKey(L, "theta", phi_theta.second);
  };

  // Get lua function
  lua_State* L = console.GetConsoleState();
  lua_getglobal(L, m_lua_function_name.c_str());

  // Error check lua function
  if (not lua_isfunction(L, -1))
    throw std::logic_error(fname + " attempted to access lua-function, " + m_lua_function_name +
                           ", but it seems the function"
                           " could not be retrieved.");

  // Push arguments
  LuaPush(L, cell_global_id);
  LuaPush(L, cell_material_id);

  LuaPush(L, face_node_location);
  LuaPush(L, face_node_normal);

  LuaPush(L, quadrature_angle_indices);

  {
    lua_newtable(L);
    int n = 0;
    for (auto& omega : quadrature_angle_vectors)
    {
      LuaPush(L, n + 1);
      LuaPush(L, omega);
      lua_settable(L, -3);
      ++n;
    }
  } // push omegas

  {
    lua_newtable(L);
    int n = 0;
    for (auto& phi_theta : quadrature_phi_theta_angles)
    {
      LuaPush(L, n + 1);
      PushPhiThetaPairTable(L, phi_theta);
      lua_settable(L, -3);
      ++n;
    }
  } // push phi_theta_pairs

  LuaPush(L, group_indices);

  LuaPush(L, time);

  std::vector<double> psi;
  // 9 arguments, 1 result (table), 0=original error object
  if (lua_pcall(L, 9, 1, 0) == 0)
  {
    LuaCheckTableValue(fname, L, -1);
    size_t table_length = lua_rawlen(L, -1);
    psi.reserve(table_length);
    for (size_t i = 0; i < table_length; ++i)
    {
      LuaPush(L, static_cast<lua_Integer>(i) + 1);
      lua_gettable(L, -2);
      psi.push_back(lua_tonumber(L, -1));
      lua_pop(L, 1);
    }
  }
  else
    throw std::logic_error(fname + " attempted to call lua-function, " + m_lua_function_name +
                           ", but the call failed.");

  lua_pop(L, 1); // pop the table, or error code

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
