// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include "framework/graphs/graph_partitioner.h"
#include "framework/parameters/input_parameters.h"

namespace opensn
{

class PETScGraphPartitioner : public GraphPartitioner
{
public:
  explicit PETScGraphPartitioner(const InputParameters& params);

  std::vector<int64_t> Partition(const std::vector<std::vector<uint64_t>>& graph,
                                 const std::vector<Vector3>& centroids,
                                 int number_of_parts) override;

protected:
  const std::string type_;

public:
  static InputParameters GetInputParameters();
  static std::shared_ptr<PETScGraphPartitioner> Create(const ParameterBlock& params);
};

} // namespace opensn
