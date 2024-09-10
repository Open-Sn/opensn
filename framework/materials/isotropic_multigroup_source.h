// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include "framework/materials/material_property.h"
#include <memory>

namespace opensn
{

/// Basic isotropic multi-group source material property.
class IsotropicMultiGroupSource : public MaterialProperty
{
public:
  std::vector<double> source_value_g;

  IsotropicMultiGroupSource() : MaterialProperty(PropertyType::ISOTROPIC_MG_SOURCE) {}

public:
  static std::shared_ptr<IsotropicMultiGroupSource> FromArray(const std::vector<double>& values)
  {
    auto mg_src = std::make_shared<IsotropicMultiGroupSource>();
    mg_src->source_value_g = values;
    return mg_src;
  }
};

} // namespace opensn
