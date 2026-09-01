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

/**
 * Mesh generator that extrudes 2D geometry into 3D by stacking layers.
 *
 * Extrusion layers are specified using an ExtrusionLayer specification which takes either pairs of
 * parameters: Pair A = "n" and "z", or Pair B = "n" and "h". When pair A is used then the z-levels
 * will be computed automatically. Vice versa, when pair B is used then the h-levels will be
 * computed automatically. Layers can be specified with a mixture of Pair A and Pair B. For example:
 * Two main layers, one specified using a height, and the other specified using a z-level.
 */
class ExtruderMeshGenerator : public MeshGenerator
{
public:
  explicit ExtruderMeshGenerator(const InputParameters& params);

  std::shared_ptr<UnpartitionedMesh>
  GenerateUnpartitionedMesh(std::shared_ptr<UnpartitionedMesh> input_umesh) override;

protected:
  const std::string top_boundary_name_;
  const std::string bottom_boundary_name_;

  std::vector<ExtrusionLayer> layers_;

public:
  static InputParameters GetInputParameters();
  static std::shared_ptr<ExtruderMeshGenerator> Create(const ParameterBlock& params);
};

} // namespace opensn
