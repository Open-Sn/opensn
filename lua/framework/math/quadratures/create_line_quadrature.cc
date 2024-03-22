#include "framework/lua.h"

#include "framework/runtime.h"
#include "framework/logging/log.h"

#include "framework/object_factory.h"

#include "quadratures_lua.h"
#include "framework/console/console.h"

using namespace opensn;

namespace opensnlua
{

RegisterLuaFunctionNamespace(CreateLineQuadrature, squad, CreateLineQuadrature);

int
CreateLineQuadrature(lua_State* L)
{
  const std::string fname = "squad.CreateLineQuadrature";
  LuaCheckArgs<int, int>(L, fname);

  // Parse argument
  auto ident = LuaArg<int>(L, 1);
  auto N = LuaArg<int>(L, 2);
  auto verbose = LuaArgOptional<bool>(L, 3, false);

  ParameterBlock params;
  params.AddParameter("verbose", verbose);
  params.AddParameter("N", N);

  auto& obj_factory = ObjectFactory::GetInstance();

  if (ident == 1) // GAUSS_LEGENDRE
  {
    opensn::log.Log() << "Creating Gauss-Legendre Quadrature\n";

    const size_t handle =
      obj_factory.MakeRegisteredObjectOfType("squad::QuadratureGaussLegendre", params);

    LuaPush(L, handle);
    return 1;
  }
  else if (ident == 2) // GAUSS_CHEBYSHEV
  {
    opensn::log.Log() << "Creating Gauss-Chebyshev Quadrature\n";

    const size_t handle =
      obj_factory.MakeRegisteredObjectOfType("squad::QuadratureGaussChebyshev", params);

    LuaPush(L, handle);
    return 1;
  }
  return 0;
}

} // namespace opensnlua
