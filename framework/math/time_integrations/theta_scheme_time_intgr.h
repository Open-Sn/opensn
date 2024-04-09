// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include "framework/math/time_integrations/time_integration.h"

namespace opensn
{

class ThetaSchemeTimeIntegration : public TimeIntegration
{
private:
  const double theta_;

public:
  static InputParameters GetInputParameters();
  explicit ThetaSchemeTimeIntegration(const InputParameters& params);

  double ThetaFactor() const;
};

} // namespace opensn
