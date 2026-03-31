#include "framework/math/functions/piecewise_linear_1d.h"
#include "framework/parameters/input_parameters.h"
#include <gtest/gtest.h>

using namespace opensn;

namespace
{

PiecewiseLinear1D
MakePiecewiseLinear(const std::vector<double>& x_values, const std::vector<double>& y_values)
{
  ParameterBlock params_block("params");
  params_block.AddParameter(ParameterBlock("x_values", x_values));
  params_block.AddParameter(ParameterBlock("y_values", y_values));

  auto params = PiecewiseLinear1D::GetInputParameters();
  params.AssignParameters(params_block);
  return PiecewiseLinear1D(params);
}

} // namespace

TEST(PiecewiseLinear1DTest, EvaluateAndSlope)
{
  auto func = MakePiecewiseLinear({0.0, 1.0, 3.0}, {0.0, 2.0, 2.0});

  EXPECT_NEAR(func.GetScalarFunction1Parameter(-1.0), 0.0, 1e-12);
  EXPECT_NEAR(func.GetScalarFunction1Parameter(0.5), 1.0, 1e-12);
  EXPECT_NEAR(func.GetScalarFunction1Parameter(2.0), 2.0, 1e-12);
  EXPECT_NEAR(func.GetScalarFunction1Parameter(5.0), 2.0, 1e-12);

  EXPECT_NEAR(func.GetScalarFunctionSlope1Parameter(-1.0), 0.0, 1e-12);
  EXPECT_NEAR(func.GetScalarFunctionSlope1Parameter(0.25), 2.0, 1e-12);
  EXPECT_NEAR(func.GetScalarFunctionSlope1Parameter(2.5), 0.0, 1e-12);
  EXPECT_NEAR(func.GetScalarFunctionSlope1Parameter(5.0), 0.0, 1e-12);

  auto v = func.Evaluate({0.5});
  ASSERT_EQ(v.size(), 1);
  EXPECT_NEAR(v[0], 1.0, 1e-12);

  auto s = func.EvaluateSlope({0.25});
  ASSERT_EQ(s.size(), 1);
  EXPECT_NEAR(s[0], 2.0, 1e-12);
}

TEST(PiecewiseLinear1DTest, InvalidInputs)
{
  EXPECT_THROW(MakePiecewiseLinear({0.0}, {1.0}), std::invalid_argument);
  EXPECT_THROW(MakePiecewiseLinear({0.0, 1.0}, {1.0}), std::invalid_argument);
  EXPECT_THROW(MakePiecewiseLinear({0.0, 1.0, 0.5}, {0.0, 1.0, 2.0}), std::invalid_argument);

  auto func = MakePiecewiseLinear({0.0, 1.0}, {0.0, 1.0});
  EXPECT_THROW(func.Evaluate({}), std::invalid_argument);
  EXPECT_THROW(func.Evaluate({0.0, 1.0}), std::invalid_argument);
  EXPECT_THROW(func.EvaluateSlope({}), std::invalid_argument);
  EXPECT_THROW(func.EvaluateSlope({0.0, 1.0}), std::invalid_argument);

  EXPECT_THROW(func.GetScalarFunction4Parameters(0.0, 0.0, 0.0, 0.0), std::logic_error);
}
