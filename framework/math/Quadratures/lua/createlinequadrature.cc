#include "chi_lua.h"

#include "chi_runtime.h"
#include "chi_log.h"

#include "ChiObjectFactory.h"

#include "quadratures_lua.h"
#include "console/chi_console.h"

RegisterLuaFunctionAsIs(chiCreateLineQuadrature);

// ########################################################## Create empty
// system
/** Creates a quadrature.
 *
\param QuadratureType int Quadrature identifier.
\param NumberOfPoints int Number of quadrature points.
\param VerboseFlag bool As the name implies. Default: false.

##_

###QuadratureType:\n
 GAUSS_LEGENDRE = Gauss-Legendre quadrature.
 GAUSS_CHEBYSHEV = Gauss-Chebyshev quadrature.

\return Returns a unique handle to the created quadrature rule

\ingroup LuaQuadrature
\author Jan*/
int chiCreateLineQuadrature(lua_State* L)
{
  const std::string fname = __FUNCTION__;
  const int num_args = lua_gettop(L);

  if (not((num_args == 2) or (num_args == 3)))
    LuaPostArgAmountError(fname, 2, num_args);

  LuaCheckNilValue(fname, L, 1);
  LuaCheckNilValue(fname, L, 2);

  //============================================= Parse argument
  int ident = lua_tonumber(L, 1);
  int N = lua_tonumber(L, 2);
  bool verbose = false;
  if (num_args == 3) verbose = lua_toboolean(L, 3);

  chi::ParameterBlock params;
  params.AddParameter("verbose", verbose);
  params.AddParameter("N", N);

  auto& obj_factory = ChiObjectFactory::GetInstance();

  if (ident == 1) // GAUSS_LEGENDRE
  {
    Chi::log.Log() << "Creating Gauss-Legendre Quadrature\n";

    const size_t handle = obj_factory.MakeRegisteredObjectOfType(
      "chi_math::QuadratureGaussLegendre", params);

    lua_pushinteger(L, static_cast<lua_Integer>(handle));
    return 1;
  }
  else if (ident == 2) // GAUSS_CHEBYSHEV
  {
    Chi::log.Log() << "Creating Gauss-Chebyshev Quadrature\n";

    const size_t handle = obj_factory.MakeRegisteredObjectOfType(
      "chi_math::QuadratureGaussChebyshev", params);

    lua_pushinteger(L, static_cast<lua_Integer>(handle));
    return 1;
  }
  return 0;
}