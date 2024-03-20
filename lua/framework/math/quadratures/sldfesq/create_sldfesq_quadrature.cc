#include "framework/lua.h"

#include "framework/runtime.h"

#include "framework/math/quadratures/sldfesq/sldfe_sq.h"

#include "framework/console/console.h"
#include "sldfe_lua.h"

using namespace opensn;

RegisterLuaFunctionNamespace(CreateSLDFESQAngularQuadrature, aquad, CreateSLDFESQAngularQuadrature);

int
CreateSLDFESQAngularQuadrature(lua_State* L)
{
  int num_args = lua_gettop(L);
  if (num_args != 1)
    LuaPostArgAmountError("CreateSLDFESQAngularQuadrature", 1, num_args);

  int init_refinement_level = lua_tonumber(L, 1);

  auto sldfesq = new SimplifiedLDFESQ::Quadrature;
  sldfesq->GenerateInitialRefinement(init_refinement_level);

  std::shared_ptr<AngularQuadrature> new_ang_quad =
    std::shared_ptr<SimplifiedLDFESQ::Quadrature>(sldfesq);

  opensn::angular_quadrature_stack.push_back(new_ang_quad);
  const size_t index = opensn::angular_quadrature_stack.size() - 1;
  lua_pushnumber(L, static_cast<lua_Number>(index));

  return 1;
}
