#include "chi_lua.h"

#include "chi_runtime.h"

#include "math/Quadratures/angular_product_quadrature.h"

#include "chi_log.h"

#include "quadratures_lua.h"
#include "console/chi_console.h"

RegisterLuaFunctionAsIs(chiGetProductQuadrature);

//########################################################## Get product quadrature
/** Get the values of a product quadrature

\param QuadHandle int Handle to an existing product quadrature.

\return Table A lua table with each entry being another table with entries
       .weight .polar .azimuthal.

\ingroup LuaQuadrature
\author Jan*/
int chiGetProductQuadrature(lua_State *L)
{
  int num_args = lua_gettop(L);
  if (num_args != 1)
    LuaPostArgAmountError("chiGetProductQuadrature",1,num_args);

  int handle = lua_tonumber(L,1);

  std::shared_ptr<chi_math::ProductQuadrature> quad;
  try{
    auto ang_quad = Chi::angular_quadrature_stack.at(handle);
    if (ang_quad->type_ == chi_math::AngularQuadratureType::ProductQuadrature)
      quad = std::static_pointer_cast<chi_math::ProductQuadrature>(ang_quad);
    else
    {
      Chi::log.LogAllError()
        << "chiGetProductQuadrature: Provided quadrature handle points to "
           "a quadrature that is not a product quadrature.";
      Chi::Exit(EXIT_FAILURE);
    }
  }
  catch (const std::out_of_range& o){
    Chi::log.LogAllError()
      << "chiGetProductQuadrature: Invalid quadrature handle.";
    Chi::Exit(EXIT_FAILURE);
  }

  lua_newtable(L);
  for (size_t n=0; n<quad->weights_.size(); ++n)
  {
    lua_pushnumber(L,n+1);
    lua_newtable(L);

    lua_pushstring(L,"weight");
    lua_pushnumber(L,quad->weights_[n]);
    lua_settable(L,-3);

    lua_pushstring(L,"polar");
    lua_pushnumber(L,quad->abscissae_[n].theta);
    lua_settable(L,-3);

    lua_pushstring(L,"azimuthal");
    lua_pushnumber(L,quad->abscissae_[n].phi);
    lua_settable(L,-3);

    lua_settable(L,-3);
  }

  return 1;
}