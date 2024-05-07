// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include "framework/materials/material_property.h"

namespace opensn
{

/** Basic isotropic multi-group source material property. */
class IsotropicMultiGroupSource : public MaterialProperty
{
public:
  std::vector<double> source_value_g;

  IsotropicMultiGroupSource() : MaterialProperty(PropertyType::ISOTROPIC_MG_SOURCE) {}
};

} // namespace opensn
