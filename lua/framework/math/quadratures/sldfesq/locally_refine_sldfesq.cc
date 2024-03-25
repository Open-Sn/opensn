#include "framework/lua.h"
#include "framework/math/quadratures/angular/sldfe_sq.h"
#include "framework/console/console.h"
#include "framework/logging/log.h"
#include "framework/runtime.h"
#include "sldfe_lua.h"

using namespace opensn;

namespace opensnlua
{

RegisterLuaFunctionNamespace(LocallyRefineSLDFESQAngularQuadrature, aquad, LocallyRefineSLDFESQ);

int
LocallyRefineSLDFESQAngularQuadrature(lua_State* L)
{
  int num_args = lua_gettop(L);
  if ((num_args != 3) and (num_args != 4))
    LuaPostArgAmountError("LocallyRefineSLDFESQAngularQuadrature", 3, num_args);

  int handle = lua_tonumber(L, 1);

  Vector3 ref_dir;
  if (lua_istable(L, 2))
  {
    lua_pushnumber(L, 1);
    lua_gettable(L, 2);
    ref_dir.x = lua_tonumber(L, -1);
    lua_pop(L, 1);

    lua_pushnumber(L, 2);
    lua_gettable(L, 2);
    ref_dir.y = lua_tonumber(L, -1);
    lua_pop(L, 1);

    lua_pushnumber(L, 3);
    lua_gettable(L, 2);
    ref_dir.z = lua_tonumber(L, -1);
    lua_pop(L, 1);
  }
  else
  {
    opensn::log.LogAllError() << "LocallyRefineSLDFESQAngularQuadrature: "
                                 "Second argument expected to be table {a,b,c}.";
    opensn::Exit(EXIT_FAILURE);
  }

  double cone_size = lua_tonumber(L, 3);

  bool ref_dir_as_plane_normal = false;
  if (num_args == 4)
    ref_dir_as_plane_normal = lua_toboolean(L, 4);

  try
  {
    auto ref_quadrature = opensn::angular_quadrature_stack.at(handle);
    if (ref_quadrature->type_ == AngularQuadratureType::SLDFESQ)
    {
      auto sldfesq = std::dynamic_pointer_cast<SimplifiedLDFESQ::Quadrature>(ref_quadrature);

      sldfesq->LocallyRefine(ref_dir, cone_size, ref_dir_as_plane_normal);
    }
    else
    {
      opensn::log.LogAllError() << "LocallyRefineSLDFESQAngularQuadrature: "
                                   "Invalid angular quadrature type.";
      opensn::Exit(EXIT_FAILURE);
    }
  }
  catch (const std::out_of_range& o)
  {
    opensn::log.LogAllError() << "LocallyRefineSLDFESQAngularQuadrature: "
                                 "Invalid handle to angular quadrature.";
    opensn::Exit(EXIT_FAILURE);
  }
  catch (...)
  {
    opensn::log.LogAllError() << "LocallyRefineSLDFESQAngularQuadrature: "
                                 "Call failed with unknown error.";
    opensn::Exit(EXIT_FAILURE);
  }

  return 0;
}

} // namespace opensnlua
