// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

namespace opensnlua
{

/**
 * Sets boundary numbers on boundaries orthogonal to the cardinal directions as xmax=0, xmin=1,
 * ymax=2, ymin=3, zmax=4, zmin=5.
 */
int MeshSetupOrthogonalBoundaries(lua_State* L);

} // namespace opensnlua
