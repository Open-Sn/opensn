// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include "framework/mesh/mesh_generator/mesh_generator.h"

namespace opensn
{

struct ExtrusionLayer
{
  static InputParameters GetInputParameters();

  const double height;
  const uint32_t num_sub_layers;
};

class ExtruderMeshGenerator : public MeshGenerator
{
public:
  static InputParameters GetInputParameters();
  explicit ExtruderMeshGenerator(const InputParameters& params);

protected:
  std::shared_ptr<UnpartitionedMesh>
  GenerateUnpartitionedMesh(std::shared_ptr<UnpartitionedMesh> input_umesh) override;

  const std::string top_boundary_name_;
  const std::string bottom_boundary_name_;

  std::vector<ExtrusionLayer> layers_;
};

} // namespace opensn
