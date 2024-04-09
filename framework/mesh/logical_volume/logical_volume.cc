// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "framework/mesh/logical_volume/logical_volume.h"

namespace opensn
{

InputParameters
LogicalVolume::GetInputParameters()
{
  return Object::GetInputParameters();
}

LogicalVolume::LogicalVolume(const InputParameters& params) : Object(params)
{
}

} // namespace opensn
