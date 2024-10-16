// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include "framework/materials/material_property.h"

namespace opensn
{

/// Simple scalar-valued material property.
class ScalarValue : public MaterialProperty
{
private:
  double value_ = 1.0;

public:
  ScalarValue() : MaterialProperty(PropertyType::SCALAR_VALUE) {}

  double GetScalarValue() override { return value_; }

  void Set(double value) { value_ = value; }
};

} // namespace opensn
