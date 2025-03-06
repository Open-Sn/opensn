// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include "framework/mesh/mesh_generator/mesh_generator.h"

namespace opensn
{

class OrthogonalMeshGenerator : public MeshGenerator
{
public:
  explicit OrthogonalMeshGenerator(const InputParameters& params);

protected:
  std::shared_ptr<UnpartitionedMesh>
  GenerateUnpartitionedMesh(std::shared_ptr<UnpartitionedMesh> input_umesh) override;

  std::vector<std::vector<double>> node_sets_;

public:
  static InputParameters GetInputParameters();
  static std::shared_ptr<OrthogonalMeshGenerator> Create(const ParameterBlock& params);

protected:
  static std::shared_ptr<UnpartitionedMesh>
  CreateUnpartitioned1DOrthoMesh(const std::vector<double>& vertices);

  static std::shared_ptr<UnpartitionedMesh>
  CreateUnpartitioned2DOrthoMesh(const std::vector<double>& vertices_1d_x,
                                 const std::vector<double>& vertices_1d_y);

  static std::shared_ptr<UnpartitionedMesh>
  CreateUnpartitioned3DOrthoMesh(const std::vector<double>& vertices_1d_x,
                                 const std::vector<double>& vertices_1d_y,
                                 const std::vector<double>& vertices_1d_z);
};

} // namespace opensn
