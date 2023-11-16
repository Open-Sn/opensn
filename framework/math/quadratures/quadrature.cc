#include "framework/math/quadratures/quadrature.h"

#include "framework/logging/log_exceptions.h"

#include "framework/console/console.h"

namespace opensn
{

InputParameters GetSyntax_Get1DQuadratureData();
ParameterBlock Get1DQuadratureData(const InputParameters& params);

RegisterWrapperFunction(chi_math,
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

  auto& quad = Chi::GetStackItem<Quadrature>(Chi::object_stack, handle, __FUNCTION__);

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

InputParameters
Quadrature::GetInputParameters()
{
  InputParameters params = Object::GetInputParameters();

  params.SetGeneralDescription("\\defgroup chi_math__Quadrature\n"
                               "\\ingroup LuaQuadrature\n"
                               "Base class for 1D quadratures");

  params.AddRequiredParameter<int>("order", "Quadrature order.");
  params.AddOptionalParameter("verbose", false, "Enables verbose operations");

  params.ConstrainParameterRange("order", AllowableRangeLowHighLimit::New(0, 43));

  return params;
}

Quadrature::Quadrature(const InputParameters& params)
  : Object(params),
    order_(static_cast<QuadratureOrder>(params.GetParamValue<int>("order"))),
    verbose_(params.GetParamValue<bool>("verbose"))
{
}

void
Quadrature::SetRange(const std::pair<double, double>& in_range)
{
  const auto& old_range = range_;
  const auto& new_range = in_range;

  const double h_new = new_range.second - new_range.first;
  const double h_old = old_range.second - old_range.first;

  ChiInvalidArgumentIf(h_new <= 0.0 or h_old <= 0.0, "Called with negative or zero ranges.");

  ChiInvalidArgumentIf(qpoints_.empty(), "Called with no abscissae initialized.");

  const double scale_factor = h_new / h_old;

  for (unsigned int i = 0; i < qpoints_.size(); ++i)
  {
    qpoints_[i](0) = new_range.first + (qpoints_[i][0] - old_range.first) * scale_factor;

    weights_[i] *= scale_factor;
  }

  range_ = in_range;
}

} // namespace opensn
