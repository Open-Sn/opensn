// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include "framework/materials/material_property.h"

namespace opensn
{

/** Basic thermal conductivity material property.*/
class IsotropicMultiGrpSource : public MaterialProperty
{
public:
  std::vector<double> source_value_g;

  IsotropicMultiGrpSource() : MaterialProperty(PropertyType::ISOTROPIC_MG_SOURCE) {}
};

} // namespace opensn
