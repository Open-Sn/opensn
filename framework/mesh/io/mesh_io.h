// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include "framework/mesh/unpartitioned_mesh/unpartitioned_mesh.h"
#include "framework/mesh/mesh_continuum/grid_vtk_utils.h"

namespace opensn
{

class MeshIO
{
public:
  static std::shared_ptr<UnpartitionedMesh> FromExodus(const UnpartitionedMesh::Options& options);
  static std::shared_ptr<UnpartitionedMesh> FromVTU(const UnpartitionedMesh::Options& options);
  static std::shared_ptr<UnpartitionedMesh> FromPVTU(const UnpartitionedMesh::Options& options);
  static std::shared_ptr<UnpartitionedMesh>
  FromEnsightGold(const UnpartitionedMesh::Options& options);
  static std::shared_ptr<UnpartitionedMesh> FromOBJ(const UnpartitionedMesh::Options& options);
  static std::shared_ptr<UnpartitionedMesh> FromGmsh(const UnpartitionedMesh::Options& options);
};

} // namespace opensn
