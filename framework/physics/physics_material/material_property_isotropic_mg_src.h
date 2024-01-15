#pragma once

#include "framework/physics/physics_material/material_property_base.h"

namespace opensn
{

/** Basic thermal conductivity material property.*/
class IsotropicMultiGrpSource : public MaterialProperty
{
public:
  std::vector<double> source_value_g_;

  IsotropicMultiGrpSource() : MaterialProperty(PropertyType::ISOTROPIC_MG_SOURCE) {}
};

} // namespace opensn
