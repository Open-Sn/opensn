#include "framework/chi_lua.h"

#include "framework/chi_runtime.h"

#include "framework/math/Quadratures/SLDFESQ/sldfe_sq.h"

#include "framework/console/chi_console.h"
#include "sldfe_lua.h"

RegisterLuaFunctionAsIs(chiCreateSLDFESQAngularQuadrature);

int
chiCreateSLDFESQAngularQuadrature(lua_State* L)
{
  int num_args = lua_gettop(L);
  if (num_args != 1) LuaPostArgAmountError("chiCreateSLDFESQAngularQuadrature", 1, num_args);

  int init_refinement_level = lua_tonumber(L, 1);

  auto sldfesq = new chi_math::SimplifiedLDFESQ::Quadrature;
  sldfesq->GenerateInitialRefinement(init_refinement_level);

  std::shared_ptr<chi_math::AngularQuadrature> new_ang_quad =
    std::shared_ptr<chi_math::SimplifiedLDFESQ::Quadrature>(sldfesq);

  Chi::angular_quadrature_stack.push_back(new_ang_quad);
  const size_t index = Chi::angular_quadrature_stack.size() - 1;
  lua_pushnumber(L, static_cast<lua_Number>(index));

  return 1;
}
