// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include "framework/materials/material_property.h"

#include <vector>
#include <memory>

namespace opensn
{

/** Base class for materials used in physics simulations.*/
class Material
{
public:
  std::vector<std::shared_ptr<MaterialProperty>> properties_{};
  std::string name_ = "Unnamed Material";
};

} // namespace opensn
