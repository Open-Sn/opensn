// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

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

  SteppingMethod GetMethod() const;

  virtual ~TimeIntegration() = default;
};

} // namespace opensn
