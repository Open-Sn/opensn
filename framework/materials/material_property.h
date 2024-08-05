// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include <string>
#include <vector>

namespace opensn
{

enum class PropertyType
{
  SCALAR_VALUE = 1,
  TRANSPORT_XSECTIONS = 10,
  ISOTROPIC_MG_SOURCE = 11
};

/// Base class for material properties.
class MaterialProperty
{
private:
  const PropertyType type_;

public:
  std::string property_name;

  explicit MaterialProperty(PropertyType type) : type_(type) {}

  virtual ~MaterialProperty() = default;

  PropertyType Type() { return type_; }

  virtual double GetScalarValue() { return 0.0; }
};

} // namespace opensn
