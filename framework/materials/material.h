// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include "framework/materials/isotropic_multigroup_source.h"
#include "framework/materials/material_property.h"
#include "framework/materials/multi_group_xs/multi_group_xs.h"
#include "framework/materials/isotropic_multigroup_source.h"
#include "framework/logging/log_exceptions.h"
#include <vector>
#include <memory>

namespace opensn
{

enum class OperationType
{
  SINGLE_VALUE = 0,
  FROM_ARRAY = 1,
  SIMPLE_ONE_GROUP = 20,
  EXISTING = 22,
  OPENSN_XSFILE = 23,
  OPENMC_XSLIB = 24
};

/// Base class for materials used in physics simulations.
class Material : public std::enable_shared_from_this<Material>
{
public:
  std::vector<std::shared_ptr<MaterialProperty>> properties{};
  std::string name = "Unnamed Material";

  void SetTransportXSections(std::shared_ptr<MultiGroupXS> xs);
  void SetIsotropicMGSource(std::shared_ptr<IsotropicMultiGroupSource> mg_src);

private:
  int GetPropertyIndex(PropertyType property_type);

  template <typename TYPE>
  std::shared_ptr<MaterialProperty>& GetProperty(int property_index)
  {
    // Create the property, if no location was found
    if (property_index == -1)
    {
      auto property = std::make_shared<TYPE>();
      properties.push_back(property);
      property_index = static_cast<int>(properties.size()) - 1;
    }
    OpenSnLogicalErrorIf(property_index < 0, "Error creating or finding MultiGroupXS property.");

    // Get the property
    return properties.at(property_index);
  }
};

} // namespace opensn
