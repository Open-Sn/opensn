// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include "framework/graphs/graph_partitioner.h"

namespace opensn
{

/**
 * Basic linear partitioner intended for testing only.
 *
 * Partitions cells by their linear `global_id` index rather than by graph connectivity, so
 * results are generally poor for unstructured grids.
 */
class LinearGraphPartitioner : public GraphPartitioner
{
public:
  explicit LinearGraphPartitioner(const InputParameters& params);

  std::vector<int> Partition(const std::vector<std::vector<uint64_t>>& graph,
                             const std::vector<Vector3>& centroids,
                             int number_of_parts) override;

protected:
  const int all_to_rank_;

public:
  static InputParameters GetInputParameters();
  static std::shared_ptr<LinearGraphPartitioner> Create(const ParameterBlock& params);
};

} // namespace opensn
