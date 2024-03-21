#include "framework/lua.h"
#include "framework/math/quadratures/angular/angular_quadrature.h"
#include "framework/logging/log.h"
#include "framework/console/console.h"
#include "framework/runtime.h"
#include "quadratures_lua.h"

namespace opensnlua
{

RegisterLuaFunctionNamespace(CreateCustomAngularQuadrature, aquad, CreateCustomAngularQuadrature);

int
CreateCustomAngularQuadrature(lua_State* L)
{
  const std::string fname = __FUNCTION__;
  const int num_args = lua_gettop(L);

  if (num_args != 3)
    LuaPostArgAmountError(fname, 3, num_args);

  std::vector<double> azi_angles = LuaArgVector<double>(L, 1);
  std::vector<double> pol_angles = LuaArgVector<double>(L, 2);
  std::vector<double> weights = LuaArgVector<double>(L, 3);

  if ((azi_angles.size() != pol_angles.size()) or (azi_angles.size() != weights.size()))
  {
    opensn::log.LogAllError() << fname + ": Tables lengths supplied "
                                         "are not of equal lengths.";
    opensn::Exit(EXIT_FAILURE);
  }

  opensn::log.Log() << "Creating Custom Angular Quadrature\n";

  auto new_quad =
    std::make_shared<opensn::AngularQuadratureCustom>(azi_angles, pol_angles, weights, false);

  opensn::angular_quadrature_stack.push_back(new_quad);
  size_t index = opensn::angular_quadrature_stack.size() - 1;
  LuaPush(L, index);

  return 1;
}

} // namespace opensnlua
