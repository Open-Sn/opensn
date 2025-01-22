// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include "framework/mesh/logical_volume/logical_volume.h"

namespace opensnlua
{

void MeshSetUniformMaterialID(int mat_id);
void MeshSetMaterialIDFromLogicalVolume(std::shared_ptr<opensn::LogicalVolume> lv,
                                        int mat_id,
                                        bool sense = true);
void MeshSetBoundaryIDFromLogicalVolume(std::shared_ptr<opensn::LogicalVolume> lv,
                                        const std::string& boundary_name,
                                        bool sense = true);
void MeshSetupOrthogonalBoundaries();
void MeshSetMaterialIDFromFunction(const char* lua_fname);
void MeshSetBoundaryIDFromFunction(const char* lua_fname);
void MeshExportToPVTU(const std::string& file_name);
void MeshComputeVolumePerMaterialID();

} // namespace opensnlua
