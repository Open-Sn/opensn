// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include "framework/object.h"

namespace opensn
{

/**Base class for a material property.*/
class MaterialProperty : public Object
{
private:
  const std::string name_;

public:
  static InputParameters GetInputParameters();
  explicit MaterialProperty(const InputParameters& params);

  const std::string& TextName() const;
};

} // namespace opensn
