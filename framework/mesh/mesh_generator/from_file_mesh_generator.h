// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include "framework/mesh/mesh_generator/mesh_generator.h"

namespace opensn
{

class FromFileMeshGenerator : public MeshGenerator
{
public:
  static InputParameters GetInputParameters();
  explicit FromFileMeshGenerator(const InputParameters& params);

protected:
  std::unique_ptr<UnpartitionedMesh>
  GenerateUnpartitionedMesh(std::unique_ptr<UnpartitionedMesh> input_umesh) override;
  const std::string filename_;
  const std::string material_id_fieldname_;
  const std::string boundary_id_fieldname_;
};

} // namespace opensn
