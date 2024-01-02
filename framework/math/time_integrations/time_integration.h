#pragma once

#include "framework/object.h"
#include "framework/math/math_time_stepping.h"

namespace opensn
{

class TimeIntegration : public Object
{
private:
  SteppingMethod method_;

public:
  static InputParameters GetInputParameters();
  explicit TimeIntegration(const InputParameters& params);

  SteppingMethod Method() const;

  virtual ~TimeIntegration() = default;
};

} // namespace opensn
