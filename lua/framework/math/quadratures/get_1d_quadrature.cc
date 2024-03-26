#include "framework/math/quadratures/spatial/spatial_quadrature.h"
#include "framework/parameters/parameter_block.h"
#include "lua/framework/console/console.h"
#include "framework/runtime.h"

using namespace opensn;

namespace opensnlua
{

InputParameters GetSyntax_Get1DQuadratureData();
ParameterBlock Get1DQuadratureData(const InputParameters& params);

RegisterWrapperFunctionNamespace(math,
                                 Get1DQuadratureData,
                                 GetSyntax_Get1DQuadratureData,
                                 Get1DQuadratureData);

InputParameters
GetSyntax_Get1DQuadratureData()
{
  InputParameters params;

  params.SetGeneralDescription("Lua wrapper function for getting the data from a quadrature.");
  params.SetDocGroup("LuaQuadrature");

  params.AddRequiredParameter<size_t>("arg0", "Handle to the quadrature");

  return params;
}

ParameterBlock
Get1DQuadratureData(const InputParameters& params)
{
  ParameterBlock output;

  const size_t handle = params.GetParamValue<size_t>("arg0");

  auto& quad = opensn::GetStackItem<SpatialQuadrature>(opensn::object_stack, handle, __FUNCTION__);

  ParameterBlock qpoints_block("qpoints");
  ParameterBlock weights_block("weights");
  {
    size_t k = 0;
    for (const auto& qpointXYZ : quad.qpoints_)
    {
      qpoints_block.AddParameter(std::to_string(k++),
                                 std::vector<double>{qpointXYZ.x, qpointXYZ.y, qpointXYZ.z});
    }
    k = 0;
    for (const double w : quad.weights_)
      weights_block.AddParameter(std::to_string(k++), w);
  }
  qpoints_block.ChangeToArray();
  weights_block.ChangeToArray();

  output.AddParameter(qpoints_block);
  output.AddParameter(weights_block);

  ParameterBlock output_as_table;
  output_as_table.AddParameter(output);

  return output_as_table;
}

} // namespace opensnlua
