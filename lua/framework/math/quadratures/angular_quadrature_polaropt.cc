#include "framework/chi_lua.h"

#include "framework/math/Quadratures/angular_quadrature_base.h"

#include "framework/chi_runtime.h"
#include "framework/logging/chi_log.h"

#include "quadratures_lua.h"
#include "framework/console/chi_console.h"

RegisterLuaFunctionAsIs(chiOptimizeAngularQuadratureForPolarSymmetry);

int
chiOptimizeAngularQuadratureForPolarSymmetry(lua_State* L)
{
  const std::string fname = __FUNCTION__;
  const int num_args = lua_gettop(L);

  if (num_args < 1) LuaPostArgAmountError(fname, 1, num_args);

  const int handle = lua_tointeger(L, 1);
  double normalization = -1.0;
  if (num_args == 2) normalization = lua_tonumber(L, 2);

  auto& quadrature =
    Chi::GetStackItem<chi_math::AngularQuadrature>(Chi::angular_quadrature_stack, handle, fname);

  if (normalization > 0.0)
    Chi::log.Log() << "Optimizing angular quadrature for polar symmetry. using "
                   << "normalization factor " << normalization << ".";

  quadrature.OptimizeForPolarSymmetry(normalization);

  return 0;
}
