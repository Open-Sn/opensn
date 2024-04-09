// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include "framework/object.h"

namespace opensn
{

/**A generic material object used to group together multiple properties.*/
class Material : public Object
{
private:
  std::string name_;

public:
  static InputParameters GetInputParameters();
  explicit Material(const InputParameters& params);
};

} // namespace opensn
