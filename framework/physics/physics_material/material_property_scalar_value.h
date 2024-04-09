// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include "framework/physics/physics_material/material_property_base.h"

namespace opensn
{

/**Simple scalar material property.*/
class ScalarValue : public PhysicsMaterialProperty
{
public:
  double value_ = 1.0;

  ScalarValue() : PhysicsMaterialProperty(PropertyType::SCALAR_VALUE) {}

  double GetScalarValue() override { return value_; }
};

} // namespace opensn
