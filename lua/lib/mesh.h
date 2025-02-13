// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include "framework/mesh/logical_volume/logical_volume.h"
#include "framework/mesh/mesh.h"

namespace opensnlua
{

void MeshSetMaterialIDFromFunction(std::shared_ptr<opensn::MeshContinuum> grid,
                                   const char* lua_fname);
void MeshSetBoundaryIDFromFunction(std::shared_ptr<opensn::MeshContinuum> grid,
                                   const char* lua_fname);
void MeshExportToPVTU(std::shared_ptr<opensn::MeshContinuum> grid, const std::string& file_name);

} // namespace opensnlua
