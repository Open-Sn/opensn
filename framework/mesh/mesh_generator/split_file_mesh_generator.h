#pragma once

#include "framework/mesh/mesh_generator/mesh_generator.h"

namespace opensn
{
class ByteArray;

/**Generates the mesh only on location 0, thereafter partitions the mesh
 * but instead of broadcasting the mesh to other locations it creates binary
 * mesh files for each location.*/
class SplitFileMeshGenerator : public MeshGenerator
{
public:
  static InputParameters GetInputParameters();
  explicit SplitFileMeshGenerator(const InputParameters& params);

  void Execute() override;

protected:
  void WriteSplitMesh(const std::vector<int64_t>& cell_pids,
                      const UnpartitionedMesh& umesh,
                      int num_parts);
  static void SerializeCell(const UnpartitionedMesh::LightWeightCell& cell,
                            ByteArray& serial_buffer);
  typedef std::pair<int, uint64_t> CellPIDGID;
  struct SplitMeshInfo
  {
    std::map<CellPIDGID, UnpartitionedMesh::LightWeightCell> cells_;
    std::map<uint64_t, Vector3> vertices_;
    std::map<uint64_t, std::string> boundary_id_map_;
    int mesh_attributes_;
    size_t ortho_Nx_;
    size_t ortho_Ny_;
    size_t ortho_Nz_;
    size_t num_global_vertices_;
  };
  SplitMeshInfo ReadSplitMesh();

  static std::shared_ptr<MeshContinuum> SetupLocalMesh(SplitMeshInfo& mesh_info);

  // void
  const int num_parts_;
  const std::string split_mesh_dir_path_;
  const std::string split_file_prefix_;
  const bool read_only_;
  const int verbosity_level_;
};

} // namespace opensn
