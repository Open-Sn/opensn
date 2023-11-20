#include "framework/math/quadratures/quadrature.h"

#include "framework/logging/log_exceptions.h"

namespace opensn
{

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
