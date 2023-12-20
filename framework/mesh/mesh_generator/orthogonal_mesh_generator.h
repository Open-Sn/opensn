#pragma once

#include "framework/mesh/mesh_generator/mesh_generator.h"

namespace opensn
{

class OrthogonalMeshGenerator : public MeshGenerator
{
public:
  static InputParameters GetInputParameters();
  explicit OrthogonalMeshGenerator(const InputParameters& params);

protected:
  std::unique_ptr<UnpartitionedMesh>
  GenerateUnpartitionedMesh(std::unique_ptr<UnpartitionedMesh> input_umesh) override;

  static std::unique_ptr<UnpartitionedMesh>
  CreateUnpartitioned1DOrthoMesh(const std::vector<double>& vertices);

  static std::unique_ptr<UnpartitionedMesh>
  CreateUnpartitioned2DOrthoMesh(const std::vector<double>& vertices_1d_x,
                                 const std::vector<double>& vertices_1d_y);

  static std::unique_ptr<UnpartitionedMesh>
  CreateUnpartitioned3DOrthoMesh(const std::vector<double>& vertices_1d_x,
                                 const std::vector<double>& vertices_1d_y,
                                 const std::vector<double>& vertices_1d_z);

  std::vector<std::vector<double>> node_sets_;
};

} // namespace opensn
