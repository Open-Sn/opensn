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
  static std::shared_ptr<UnpartitionedMesh> FromExodusII(const UnpartitionedMesh::Options& options);
  static std::shared_ptr<UnpartitionedMesh> FromVTU(const UnpartitionedMesh::Options& options);
  static std::shared_ptr<UnpartitionedMesh> FromPVTU(const UnpartitionedMesh::Options& options);
  static std::shared_ptr<UnpartitionedMesh>
  FromEnsightGold(const UnpartitionedMesh::Options& options);
  static std::shared_ptr<UnpartitionedMesh> FromOBJ(const UnpartitionedMesh::Options& options);
  static std::shared_ptr<UnpartitionedMesh> FromGmsh(const UnpartitionedMesh::Options& options);
  static std::shared_ptr<UnpartitionedMesh> FromOpenFOAM(const UnpartitionedMesh::Options& options);

  /**
   * Write grid cells into an OBJ file
   *
   * \param grid Grid to be stored
   * \param file_name Name of the output file
   * \param per_material Create one file per material
   */
  static void ToOBJ(const std::shared_ptr<MeshContinuum>& grid,
                    const char* file_name,
                    bool per_material = false);

  /**
   * Write grid cells into an ExodusII file
   *
   * \param grid Grid to be stored
   * \param file_name Name of the output file
   * \param per_material Create one file per material
   * \param write_node_sets Write node sets into the file
   * \param write_side_sets Write side sets into the file
   */
  static void ToExodusII(const std::shared_ptr<MeshContinuum>& grid,
                         const std::string& file_name,
                         bool write_node_sets = true,
                         bool write_side_sets = true);

  /**
   * Write grid cells into PVTU format.
   *
   * \param grid Grid to be stored
   * \param file_base_name Base name of the output file
   */
  static void ToPVTU(const std::shared_ptr<MeshContinuum>& grid, const std::string& file_base_name);

private:
  static std::shared_ptr<UnpartitionedMesh> FromGmshV41(const UnpartitionedMesh::Options& options);
  static std::shared_ptr<UnpartitionedMesh>
  FromGmshV41ASCII(const UnpartitionedMesh::Options& options);
  static std::shared_ptr<UnpartitionedMesh>
  FromGmshV41Binary(const UnpartitionedMesh::Options& options, int data_size);
  static std::shared_ptr<UnpartitionedMesh> FromGmshV22(const UnpartitionedMesh::Options& options);
};

} // namespace opensn
