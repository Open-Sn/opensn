#ifndef CHITECH_CHI_MATH_FUNCTIONS_PIECEWISELINEAR1D_H
#define CHITECH_CHI_MATH_FUNCTIONS_PIECEWISELINEAR1D_H

#include "function_dimA_to_dimB.h"

namespace chi_math::functions
{
class PiecewiseLinear1D : public FunctionDimAToDimB
{
public:
  static chi::InputParameters GetInputParameters();

  explicit PiecewiseLinear1D(const chi::InputParameters& params);

  std::vector<double> Evaluate(const std::vector<double>& values) const override;
  std::vector<double> EvaluateSlope(const std::vector<double>& values) const override;

  double ScalarFunction1Parameter(double x) const override;
  double ScalarFunctionSlope1Parameter(double x) const override;

  bool HasSlope() const override { return true; }
  bool HasCurvature() const override { return false; }

private:
  // Inputs
  /**Independent variable values.*/
  const std::vector<double> x_values_;
  /**Dependent variable values.*/
  const std::vector<double> y_values_;

  std::vector<double> slopes_;

  // Determined during construction
  /**The number of items in the discrete function values*/
  const size_t num_vals_;
  /**Distance between independent variable values. Used for interpolation.*/
  std::vector<double> delta_x_values_;
};
} // namespace chi_math::functions

#endif // CHITECH_CHI_MATH_FUNCTIONS_PIECEWISELINEAR1D_H
