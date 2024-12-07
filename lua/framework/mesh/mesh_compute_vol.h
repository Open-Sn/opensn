// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include "lua/framework/lua.h"

namespace opensnlua
{

/// Compute the volume per material ID
int MeshComputeVolumePerMaterialID(lua_State* L);

} // namespace opensnlua
