#include "framework/chi_lua.h"

#include "framework/chi_runtime.h"

#include "framework/math/Quadratures/SLDFESQ/sldfe_sq.h"

#include "framework/logging/chi_log.h"

#include "framework/mpi/chi_mpi.h"
#include "framework/console/chi_console.h"
#include "sldfe_lua.h"

RegisterLuaFunctionAsIs(chiPrintToPythonSLDFESQAngularQuadrature);

/** Outputs the quadrature information to python format.
\param handle int Handle to the reference quadrature.
\param file_name_prefix string Prefix to be used in front of file.

##_

###Example:
Example of printing a quadrature:
Example with refinement level 2 and a triple directional refinement:
\code
pquad = chiCreateSLDFESQAngularQuadrature(2)
chiLocallyRefineSLDFESQAngularQuadrature(pquad,{1,0,0},45.0*math.pi/180,false)
chiLocallyRefineSLDFESQAngularQuadrature(pquad,{1,0,0},23.0*math.pi/180,false)
chiLocallyRefineSLDFESQAngularQuadrature(pquad,{1,0,0},12.0*math.pi/180,false)
chiPrintToPythonSLDFESQAngularQuadrature(pquad,"YQuad_");
\endcode

\image html "SLDFESQr.png" width=500px

\ingroup LuaSLDFESQ
\author Jan */
int
chiPrintToPythonSLDFESQAngularQuadrature(lua_State* L)
{
  int num_args = lua_gettop(L);
  if (num_args != 2) LuaPostArgAmountError("chiPrintToPythonSLDFESQAngularQuadrature", 2, num_args);

  int handle = lua_tonumber(L, 1);
  const char* file_name = lua_tostring(L, 2);

  try
  {
    auto ref_quadrature = Chi::angular_quadrature_stack.at(handle);
    if (ref_quadrature->type_ == chi_math::AngularQuadratureType::SLDFESQ)
    {
      auto sldfesq =
        std::dynamic_pointer_cast<chi_math::SimplifiedLDFESQ::Quadrature>(ref_quadrature);

      if (Chi::mpi.location_id == 0)
      {
        sldfesq->output_filename_prefix_ = file_name;
        sldfesq->PrintQuadratureToFile();
      }
    }
    else
    {
      Chi::log.LogAllError() << "chiPrintToPythonSLDFESQAngularQuadrature: "
                                "Invalid angular quadrature type.";
      Chi::Exit(EXIT_FAILURE);
    }
  }
  catch (const std::out_of_range& o)
  {
    Chi::log.LogAllError() << "chiPrintToPythonSLDFESQAngularQuadrature: "
                              "Invalid handle to angular quadrature.";
    Chi::Exit(EXIT_FAILURE);
  }
  catch (...)
  {
    Chi::log.LogAllError() << "chiPrintToPythonSLDFESQAngularQuadrature: "
                              "Call failed with unknown error.";
    Chi::Exit(EXIT_FAILURE);
  }

  return 0;
}
