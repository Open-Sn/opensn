#include "framework/lua.h"

#include "framework/math/quadratures/angular_quadrature_base.h"

#include "framework/runtime.h"
#include "framework/logging/log.h"

#include "quadratures_lua.h"
#include "framework/console/console.h"

using namespace opensn;

namespace opensnlua
{

RegisterLuaFunctionNamespace(OptimizeAngularQuadratureForPolarSymmetry,
                             aquad,
                             OptimizeForPolarSymmetry);

int
OptimizeAngularQuadratureForPolarSymmetry(lua_State* L)
{
  const std::string fname = __FUNCTION__;
  const int num_args = lua_gettop(L);

  if (num_args < 1)
    LuaPostArgAmountError(fname, 1, num_args);

  const int handle = lua_tointeger(L, 1);
  double normalization = -1.0;
  if (num_args == 2)
    normalization = lua_tonumber(L, 2);

  auto& quadrature =
    opensn::GetStackItem<AngularQuadrature>(opensn::angular_quadrature_stack, handle, fname);

  if (normalization > 0.0)
    opensn::log.Log() << "Optimizing angular quadrature for polar symmetry. using "
                      << "normalization factor " << normalization << ".";

  quadrature.OptimizeForPolarSymmetry(normalization);

  return 0;
}

} // namespace opensnlua
