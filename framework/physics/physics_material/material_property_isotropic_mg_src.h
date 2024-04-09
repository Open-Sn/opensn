// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include "framework/physics/physics_material/material_property_base.h"

namespace opensn
{

/** Basic thermal conductivity material property.*/
class IsotropicMultiGrpSource : public PhysicsMaterialProperty
{
public:
  std::vector<double> source_value_g_;

  IsotropicMultiGrpSource() : PhysicsMaterialProperty(PropertyType::ISOTROPIC_MG_SOURCE) {}
};

} // namespace opensn
