// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include "framework/parameters/input_parameters.h"
#include <vector>

namespace opensn
{
struct Vector3;

/// Abstract base class for all partitioners
class GraphPartitioner
{
public:
  virtual ~GraphPartitioner() = default;

  /// Given a graph. Returns the partition ids of each row in the graph.
  virtual std::vector<int> Partition(const std::vector<std::vector<uint64_t>>& graph,
                                     const std::vector<Vector3>& centroids,
                                     int number_of_parts) = 0;

protected:
  explicit GraphPartitioner(const InputParameters& params);

protected:
  static InputParameters GetInputParameters();
};

} // namespace opensn
