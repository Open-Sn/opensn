#include "framework/chi_lua.h"

#include "framework/chi_runtime.h"

#include "framework/math/Quadratures/SLDFESQ/sldfe_sq.h"

#include "framework/logging/chi_log.h"
#include "framework/console/chi_console.h"
#include "sldfe_lua.h"

RegisterLuaFunctionAsIs(chiLocallyRefineSLDFESQAngularQuadrature);

// ###################################################################
/** Applies a local refinement of angles.
\param handle int. Handle to the reference quadrature.
\param reference_direction vec3 Reference vector. \f$ \vec{r} \f$
\param cone_size double Cone size in radians. \f$ \theta \f$
\param invert_logic bool Optional[Default:false]. If supplied, interprets
 SQ-splitting as when \f$|\omega \cdot \vec{r}| < \sin(\theta) \f$. Otherwise,
 SQs will be split if \f$ \omega \cdot \vec{r} > \cos(\theta)\f$

##_

###Example:
Example with refinement level 2 and a triple directional refinement:
\code
pquad = chiCreateSLDFESQAngularQuadrature(2)
chiLocallyRefineSLDFESQAngularQuadrature(pquad,{1,0,0},45.0*math.pi/180,false)
chiLocallyRefineSLDFESQAngularQuadrature(pquad,{1,0,0},23.0*math.pi/180,false)
chiLocallyRefineSLDFESQAngularQuadrature(pquad,{1,0,0},12.0*math.pi/180,false)
\endcode

\image html "SLDFESQr.png" width=500px

Example with refinement level 2 and a triple planar refinement:
\code
pquad = chiCreateSLDFESQAngularQuadrature(2)
chiLocallyRefineSLDFESQAngularQuadrature(pquad,{1,0,0},22.50*math.pi/180,true)
chiLocallyRefineSLDFESQAngularQuadrature(pquad,{1,0,0},11.75*math.pi/180,true)
chiLocallyRefineSLDFESQAngularQuadrature(pquad,{1,0,0},5.000*math.pi/180,true)
\endcode

\image html "SLDFESQp.png" width=500px

\ingroup LuaSLDFESQ
\author Jan */
int
chiLocallyRefineSLDFESQAngularQuadrature(lua_State* L)
{
  int num_args = lua_gettop(L);
  if ((num_args != 3) and (num_args != 4))
    LuaPostArgAmountError("chiLocallyRefineSLDFESQAngularQuadrature", 3, num_args);

  int handle = lua_tonumber(L, 1);

  chi_mesh::Vector3 ref_dir;
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
    Chi::log.LogAllError() << "chiLocallyRefineSLDFESQAngularQuadrature: "
                              "Second argument expected to be table {a,b,c}.";
    Chi::Exit(EXIT_FAILURE);
  }

  double cone_size = lua_tonumber(L, 3);

  bool ref_dir_as_plane_normal = false;
  if (num_args == 4) ref_dir_as_plane_normal = lua_toboolean(L, 4);

  try
  {
    auto ref_quadrature = Chi::angular_quadrature_stack.at(handle);
    if (ref_quadrature->type_ == chi_math::AngularQuadratureType::SLDFESQ)
    {
      auto sldfesq =
        std::dynamic_pointer_cast<chi_math::SimplifiedLDFESQ::Quadrature>(ref_quadrature);

      sldfesq->LocallyRefine(ref_dir, cone_size, ref_dir_as_plane_normal);
    }
    else
    {
      Chi::log.LogAllError() << "chiLocallyRefineSLDFESQAngularQuadrature: "
                                "Invalid angular quadrature type.";
      Chi::Exit(EXIT_FAILURE);
    }
  }
  catch (const std::out_of_range& o)
  {
    Chi::log.LogAllError() << "chiLocallyRefineSLDFESQAngularQuadrature: "
                              "Invalid handle to angular quadrature.";
    Chi::Exit(EXIT_FAILURE);
  }
  catch (...)
  {
    Chi::log.LogAllError() << "chiLocallyRefineSLDFESQAngularQuadrature: "
                              "Call failed with unknown error.";
    Chi::Exit(EXIT_FAILURE);
  }

  return 0;
}
