#include "framework/math/quadratures/spatial/spatial_quadrature.h"

namespace opensn
{

InputParameters
SpatialQuadrature::GetInputParameters()
{
  InputParameters params = Object::GetInputParameters();

  params.SetGeneralDescription("\\defgroup math__Quadrature\n"
                               "\\ingroup LuaQuadrature\n"
                               "Base class for spatial quadratures");

  params.AddRequiredParameter<int>("order", "Quadrature order.");

  params.AddOptionalParameter("verbose", false, "Enables verbose operations");

  return params;
}

} // namespace opensn
