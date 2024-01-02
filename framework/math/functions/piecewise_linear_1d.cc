#include "framework/math/functions/piecewise_linear_1d.h"

#include "framework/object_factory.h"

namespace opensn
{

OpenSnRegisterObject(chi_math::functions, PiecewiseLinear1D);

InputParameters
PiecewiseLinear1D::GetInputParameters()
{
  InputParameters params = FunctionDimAToDimB::GetInputParameters();

  // clang-format off
  params.SetGeneralDescription("Piecewise linear function");
  params.SetDocGroup("DocMathFunctions");
  // clang-format on

  params.AddRequiredParameterArray("x_values", "The x-values used in the interpolation function.");
  params.AddRequiredParameterArray("y_values", "The x-values used in the interpolation function.");

  params.ChangeExistingParamToOptional("input_dimension", size_t{1});
  params.ChangeExistingParamToOptional("output_dimension", size_t{1});

  return params;
}

PiecewiseLinear1D::PiecewiseLinear1D(const InputParameters& params)
  : FunctionDimAToDimB(params),
    x_values_(params.GetParamVectorValue<double>("x_values")),
    y_values_(params.GetParamVectorValue<double>("y_values")),
    num_vals_(x_values_.size())
{
  ChiInvalidArgumentIf(y_values_.size() != num_vals_,
                       std::string("Number of y-values (") + std::to_string(y_values_.size()) +
                         ") must match number of x-values (" + std::to_string(num_vals_) + ".");

  ChiInvalidArgumentIf(y_values_.size() < 2, "Number of pairs must at least be 2");

  delta_x_values_.assign(num_vals_ - 1, 0.0);
  slopes_.assign(num_vals_ - 1, 0.0);
  const size_t max_k = num_vals_ - 1;
  for (size_t k = 0; k < max_k; ++k)
  {
    delta_x_values_[k] = x_values_[k + 1] - x_values_[k];
    ChiInvalidArgumentIf(delta_x_values_[k] <= 0.0,
                         "x-values not monitonically "
                         "increasing");
    slopes_[k] = (y_values_[k + 1] - y_values_[k]) / delta_x_values_[k];
  }
}

std::vector<double>
PiecewiseLinear1D::Evaluate(const std::vector<double>& values) const
{
  if (values.size() != 1) ChiInvalidArgument("Can only be called with 1 argument.");

  return {ScalarFunction1Parameter(values.front())};
}

std::vector<double>
PiecewiseLinear1D::EvaluateSlope(const std::vector<double>& values) const
{
  if (values.size() != 1) ChiInvalidArgument("Can only be called with 1 argument.");

  return {ScalarFunctionSlope1Parameter(values.front())};
}

double
PiecewiseLinear1D::ScalarFunction1Parameter(double x) const
{
  if (x < x_values_.front()) return y_values_.front();

  if (x >= x_values_.back()) return y_values_.back();

  const size_t max_k = num_vals_ - 1;
  for (size_t k = 0; k < max_k; ++k)
    if ((x >= x_values_[k]) and (x < x_values_[k + 1]))
    {
      return (y_values_[k] * (x_values_[k + 1] - x) + y_values_[k + 1] * (x - x_values_[k])) /
             delta_x_values_[k];
    }

  ChiLogicalError("Bad trouble");
}

double
PiecewiseLinear1D::ScalarFunctionSlope1Parameter(double x) const
{
  if (x < x_values_.front()) return 0.0;

  if (x >= x_values_.back()) return 0.0;

  const size_t max_k = num_vals_ - 1;
  for (size_t k = 0; k < max_k; ++k)
    if ((x >= x_values_[k]) and (x < x_values_[k + 1])) { return slopes_[k]; }

  ChiLogicalError("Bad trouble");
}

} // namespace opensn
