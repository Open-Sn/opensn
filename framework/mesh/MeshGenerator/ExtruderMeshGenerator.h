#pragma once

#include "MeshGenerator.h"

namespace chi_mesh
{

struct ExtrusionLayer
{
  static chi::InputParameters GetInputParameters();

  const double height_;
  const uint32_t num_sub_layers_;
};

class ExtruderMeshGenerator : public MeshGenerator
{
public:
  static chi::InputParameters GetInputParameters();
  explicit ExtruderMeshGenerator(const chi::InputParameters& params);

protected:
  std::unique_ptr<UnpartitionedMesh>
  GenerateUnpartitionedMesh(std::unique_ptr<UnpartitionedMesh> input_umesh) override;

  const std::string top_boundary_name_;
  const std::string bottom_boundary_name_;

  std::vector<ExtrusionLayer> layers_;
};

} // namespace chi_mesh


