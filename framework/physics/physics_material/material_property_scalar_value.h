#pragma once

#include "framework/physics/physics_material/material_property_base.h"

namespace opensn
{

/**Simple scalar material property.*/
class ScalarValue : public MaterialProperty
{
public:
  double value_ = 1.0;

  ScalarValue() : MaterialProperty(PropertyType::SCALAR_VALUE) {}

  double GetScalarValue() override { return value_; }
};

} // namespace opensn
