// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "framework/materials/material.h"

namespace opensn
{

void
Material::SetTransportXSections(std::shared_ptr<MultiGroupXS> xs)
{
  auto property_index = GetPropertyIndex(PropertyType::TRANSPORT_XSECTIONS);
  auto& property = GetProperty<MultiGroupXS>(property_index);
  property = xs;
}

int
Material::GetPropertyIndex(PropertyType property_type)
{
  // Find the index of the material property on the material
  int property_index = -1;
  for (int p = 0; p < properties.size(); ++p)
  {
    const auto& property = properties.at(p);
    if (property->GetType() == property_type)
    {
      property_index = p;
      break;
    }
  }
  return property_index;
}

} // namespace opensn
