#include "chi_lua.h"

#include "chi_runtime.h"

#include "../sldfe_sq.h"

#include "console/chi_console.h"
#include "sldfe_lua.h"

RegisterLuaFunctionAsIs(chiCreateSLDFESQAngularQuadrature);

//###################################################################
/** Creates a Simplified Linear Discontinuous Finite Element (SLDFE)
quadrature based on Spherical Quadrilaterals (SQ). Hence SLDFE-SQ.
\param initial_refinement_level int Initial refinement level, \f$n\f$ to
       be used. The total number of angles will be \f$ 8{\times}12(n+1)^2 \f$.

##_

###Example:
Example with refinement level 2.
\code
pquad = chiCreateSLDFESQAngularQuadrature(2)
\endcode

\image html "SLDFESQBasen2.png" width=500px
With direction points:
\image html "SLDFESQBasen2Oct2.png" width=500px

\ingroup LuaSLDFESQ
\author Jan */
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