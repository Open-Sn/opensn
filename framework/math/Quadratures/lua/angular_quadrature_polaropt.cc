#include "chi_lua.h"

#include "math/Quadratures/angular_quadrature_base.h"

#include "chi_runtime.h"
#include "chi_log.h"

#include "quadratures_lua.h"
#include "console/chi_console.h"

RegisterLuaFunctionAsIs(chiOptimizeAngularQuadratureForPolarSymmetry);

//###################################################################
/**Optimizes the indicated angular quadrature for polar symmetry.
 *
\param handle        int. Handle to the quadrature to be optimized.
\param normalization double. (Optional) The normalization to be applied to the
                     modified quadrature. Any negative number will inhibit
                     renormalization. [Default=-1.0]

 ## _

 ###Example:
 Example:
\code
pquad = chiCreateProductQuadrature(GAUSS_LEGENDRE_CHEBYSHEV,2, 1)
chiOptimizeAngularQuadratureForPolarSymmetry(pqaud, 4.0*math.pi)
\endcode

\ingroup LuaQuadrature
\author Jan */
int chiOptimizeAngularQuadratureForPolarSymmetry(lua_State* L)
{
  const std::string fname = __FUNCTION__;
  const int num_args = lua_gettop(L);

  if (num_args < 1)
    LuaPostArgAmountError(fname, 1, num_args);

  const int handle = lua_tointeger(L, 1);
  double normalization = -1.0;
  if (num_args == 2)
    normalization = lua_tonumber(L, 2);

  auto& quadrature = Chi::GetStackItem<chi_math::AngularQuadrature>(
    Chi::angular_quadrature_stack, handle, fname);

  if (normalization > 0.0)
    Chi::log.Log() << "Optimizing angular quadrature for polar symmetry. using "
                   << "normalization factor " << normalization << ".";

  quadrature.OptimizeForPolarSymmetry(normalization);

  return 0;
}