// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include "framework/materials/material.h"
#include "framework/materials/multi_group_xs/multi_group_xs.h"

namespace opensnlua
{

std::shared_ptr<opensn::Material> MatAddMaterial(const std::string& name);

} // namespace opensnlua
