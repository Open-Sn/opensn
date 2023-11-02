#ifndef CHITECH_CHI_MATH_FUNCTION_DIMA_TO_DIMB_H
#define CHITECH_CHI_MATH_FUNCTION_DIMA_TO_DIMB_H

#include "ChiObject.h"
#include <functional>

namespace chi_math
{
typedef std::function<double(double)> ScalarScalarFunction;
typedef std::function<double(double, double, double, double)>
  ScalarXYZTFunction;
class FunctionDimAToDimB : public ChiObject
{
private:
  const size_t input_dimension_;
  const size_t output_dimension_;

public:
  static chi::InputParameters GetInputParameters();
  explicit FunctionDimAToDimB(const chi::InputParameters& params);

  size_t InputDimension() const { return input_dimension_; }
  size_t OutputDimension() const { return output_dimension_; }

  virtual bool HasSlope() const = 0;
  virtual bool HasCurvature() const = 0;

  virtual double ScalarFunction1Parameter(double) const;
  virtual double ScalarFunctionSlope1Parameter(double) const;
  virtual double ScalarFunctionCurvature1Parameter(double) const;

  virtual double
  ScalarFunction4Parameters(double, double, double, double) const;
  virtual double
  ScalarFunctionSlope4Parameters(double, double, double, double) const;
  virtual double
  ScalarFunctionCurvature4Parameters(double, double, double, double) const;

  virtual std::vector<double>
  Evaluate(const std::vector<double>& vals) const = 0;
  virtual std::vector<double>
  EvaluateSlope(const std::vector<double>& vals) const
  {
    return {0.0};
  }
};

} // namespace chi_math

#endif // CHITECH_CHI_MATH_FUNCTION_DIMA_TO_DIMB_H
