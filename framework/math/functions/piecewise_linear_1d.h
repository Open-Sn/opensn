// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include "framework/math/functions/function_dimA_to_dimB.h"

namespace opensn
{

class PiecewiseLinear1D : public FunctionDimAToDimB
{
public:
  static InputParameters GetInputParameters();

  explicit PiecewiseLinear1D(const InputParameters& params);

  std::vector<double> Evaluate(const std::vector<double>& values) const override;
  std::vector<double> EvaluateSlope(const std::vector<double>& values) const override;

  double GetScalarFunction1Parameter(double x) const override;
  double GetScalarFunctionSlope1Parameter(double x) const override;

  bool HasSlope() const override { return true; }
  bool HasCurvature() const override { return false; }

private:
  /// Independent variable values.
  const std::vector<double> x_values_;
  /// Dependent variable values.
  const std::vector<double> y_values_;

  std::vector<double> slopes_;

  /// The number of items in the discrete function values
  const size_t num_vals_;
  /// Distance between independent variable values. Used for interpolation.
  std::vector<double> delta_x_values_;
};

} // namespace opensn
