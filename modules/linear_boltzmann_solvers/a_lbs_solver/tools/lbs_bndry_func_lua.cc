#ifdef OPENSN_WITH_LUA
#include "modules/linear_boltzmann_solvers/a_lbs_solver/tools/lbs_bndry_func_lua.h"

#include "framework/lua.h"

#include "framework/runtime.h"
#include "framework/logging/log.h"
#include "framework/console/console.h"

namespace opensn
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
  // Utility lambdas
  auto PushVector3AsTable = [](lua_State* L, const Vector3& vec)
  {
    lua_newtable(L);

    lua_pushstring(L, "x");
    lua_pushnumber(L, vec.x);
    lua_settable(L, -3);

    lua_pushstring(L, "y");
    lua_pushnumber(L, vec.y);
    lua_settable(L, -3);

    lua_pushstring(L, "z");
    lua_pushnumber(L, vec.z);
    lua_settable(L, -3);
  };

  auto PushVecIntAsTable = [](lua_State* L, const std::vector<int>& vec)
  {
    lua_newtable(L);

    for (int i = 0; i < static_cast<int>(vec.size()); ++i)
    {
      lua_pushinteger(L, i + 1);
      lua_pushinteger(L, static_cast<lua_Integer>(vec[i]));
      lua_settable(L, -3);
    }
  };

  auto PushPhiThetaPairTable = [](lua_State* L, const std::pair<double, double>& phi_theta)
  {
    lua_newtable(L);

    lua_pushstring(L, "phi");
    lua_pushnumber(L, phi_theta.first);
    lua_settable(L, -3);

    lua_pushstring(L, "theta");
    lua_pushnumber(L, phi_theta.second);
    lua_settable(L, -3);
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
  lua_pushinteger(L, static_cast<lua_Integer>(cell_global_id));
  lua_pushinteger(L, static_cast<lua_Integer>(cell_material_id));

  PushVector3AsTable(L, face_node_location);
  PushVector3AsTable(L, face_node_normal);

  PushVecIntAsTable(L, quadrature_angle_indices);

  {
    lua_newtable(L);
    int n = 0;
    for (auto& omega : quadrature_angle_vectors)
    {
      lua_pushinteger(L, n + 1);
      PushVector3AsTable(L, omega);
      lua_settable(L, -3);
      ++n;
    }
  } // push omegas

  {
    lua_newtable(L);
    int n = 0;
    for (auto& phi_theta : quadrature_phi_theta_angles)
    {
      lua_pushinteger(L, n + 1);
      PushPhiThetaPairTable(L, phi_theta);
      lua_settable(L, -3);
      ++n;
    }
  } // push phi_theta_pairs

  PushVecIntAsTable(L, group_indices);

  lua_pushnumber(L, time);

  std::vector<double> psi;
  // 9 arguments, 1 result (table), 0=original error object
  if (lua_pcall(L, 9, 1, 0) == 0)
  {
    LuaCheckTableValue(fname, L, -1);
    size_t table_length = lua_rawlen(L, -1);
    psi.reserve(table_length);
    for (size_t i = 0; i < table_length; ++i)
    {
      lua_pushinteger(L, static_cast<lua_Integer>(i) + 1);
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
} // namespace opensn

#endif
