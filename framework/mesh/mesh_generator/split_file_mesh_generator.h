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
public:
  explicit SplitFileMeshGenerator(const InputParameters& params);

  void Execute() override;

protected:
  void WriteSplitMesh(const std::vector<int64_t>& cell_pids,
                      const UnpartitionedMesh& umesh,
                      int num_parts);
  static void SerializeCell(const UnpartitionedMesh::LightWeightCell& cell,
                            ByteArray& serial_buffer);
  struct SplitMeshInfo
  {
    unsigned int dimension;
    std::map<std::pair<int, uint64_t>, UnpartitionedMesh::LightWeightCell> cells;
    std::map<uint64_t, Vector3> vertices;
    std::map<uint64_t, std::string> boundary_id_map;
    MeshType mesh_type;
    bool extruded;
    OrthoMeshAttributes ortho_attributes;
    size_t num_global_vertices;
  };
  SplitMeshInfo ReadSplitMesh();

  static std::shared_ptr<MeshContinuum> SetupLocalMesh(SplitMeshInfo& mesh_info);

  // void
  const int num_parts_;
  const std::string split_mesh_dir_path_;
  const std::string file_prefix_;
  const bool read_only_;
  const int verbosity_level_;

public:
  static InputParameters GetInputParameters();
  static std::shared_ptr<SplitFileMeshGenerator> Create(const ParameterBlock& params);
};

} // namespace opensn
