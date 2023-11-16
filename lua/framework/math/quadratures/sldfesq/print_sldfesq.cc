#include "framework/lua.h"

#include "framework/runtime.h"

#include "framework/math/quadratures/sldfesq/sldfe_sq.h"

#include "framework/logging/log.h"

#include "framework/mpi/mpi.h"
#include "framework/console/console.h"
#include "sldfe_lua.h"

using namespace opensn;

RegisterLuaFunctionAsIs(chiPrintToPythonSLDFESQAngularQuadrature);

int
chiPrintToPythonSLDFESQAngularQuadrature(lua_State* L)
{
  int num_args = lua_gettop(L);
  if (num_args != 2) LuaPostArgAmountError("chiPrintToPythonSLDFESQAngularQuadrature", 2, num_args);

  int handle = lua_tonumber(L, 1);
  const char* file_name = lua_tostring(L, 2);

  try
  {
    auto ref_quadrature = opensn::Chi::angular_quadrature_stack.at(handle);
    if (ref_quadrature->type_ == AngularQuadratureType::SLDFESQ)
    {
      auto sldfesq = std::dynamic_pointer_cast<SimplifiedLDFESQ::Quadrature>(ref_quadrature);

      if (opensn::mpi.location_id == 0)
      {
        sldfesq->output_filename_prefix_ = file_name;
        sldfesq->PrintQuadratureToFile();
      }
    }
    else
    {
      opensn::log.LogAllError() << "chiPrintToPythonSLDFESQAngularQuadrature: "
                                   "Invalid angular quadrature type.";
      opensn::Chi::Exit(EXIT_FAILURE);
    }
  }
  catch (const std::out_of_range& o)
  {
    opensn::log.LogAllError() << "chiPrintToPythonSLDFESQAngularQuadrature: "
                                 "Invalid handle to angular quadrature.";
    opensn::Chi::Exit(EXIT_FAILURE);
  }
  catch (...)
  {
    opensn::log.LogAllError() << "chiPrintToPythonSLDFESQAngularQuadrature: "
                                 "Call failed with unknown error.";
    opensn::Chi::Exit(EXIT_FAILURE);
  }

  return 0;
}
