// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include "framework/mesh/mesh_generator/mesh_generator.h"

namespace opensn
{
class ByteArray;

/**
 * Generates the mesh only on location 0, thereafter partitions the mesh but instead of broadcasting
 * the mesh to other locations it creates binary mesh files for each location.
 */
class SplitFileMeshGenerator : public MeshGenerator
{
protected:
  struct SplitMeshInfo
  {
    unsigned int dimension;
    CoordinateSystemType coord_sys;
    MeshType mesh_type;
    bool extruded;
    OrthoMeshAttributes ortho_attributes;

    std::map<std::pair<int, uint64_t>, UnpartitionedMesh::LightWeightCell> cells;
    std::map<uint64_t, Vector3> vertices;
    std::map<uint64_t, std::string> boundary_id_map;
    size_t num_global_vertices;
  };

public:
  explicit SplitFileMeshGenerator(const InputParameters& params);

  std::shared_ptr<MeshContinuum> Execute() override;

protected:
  void WriteSplitMesh(const std::vector<int>& cell_pids,
                      const UnpartitionedMesh& umesh,
                      int num_partitions) const;

  SplitMeshInfo ReadSplitMesh() const;

  const int num_partitions_;
  const std::string split_mesh_dir_path_;
  const std::string file_prefix_;
  const bool read_only_;
  const int verbosity_level_;

public:
  static InputParameters GetInputParameters();
  static std::shared_ptr<SplitFileMeshGenerator> Create(const ParameterBlock& params);

protected:
  static std::shared_ptr<MeshContinuum> SetupLocalMesh(SplitMeshInfo& mesh_info);

  static void SerializeCell(const UnpartitionedMesh::LightWeightCell& cell,
                            ByteArray& serial_buffer);
};

} // namespace opensn
