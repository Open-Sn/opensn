// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include "framework/parameters/input_parameters.h"
#include <cstdint>
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

  /**
   * Gives every empty partition cells so that each partition has at least one.
   *
   * When no partition is empty, nothing changes. Each empty partition takes about
   * ``num_cells / num_partitions`` cells from the currently largest partition, which keeps at least
   * one cell. The cells taken are connected. They are collected in breadth-first order within the
   * donor partition, starting from the donor cell farthest from the donor's lowest-index cell.
   * Nothing changes when there are fewer cells than partitions.
   *
   * \param graph Cell adjacency, giving the neighbor indices of each cell.
   * \param cell_pids Partition id of each cell, updated in place.
   * \param num_partitions Number of partitions.
   */
  static void FillEmptyPartitions(const std::vector<std::vector<uint64_t>>& graph,
                                  std::vector<int>& cell_pids,
                                  int num_partitions);

protected:
  explicit GraphPartitioner(const InputParameters& params);

protected:
  static InputParameters GetInputParameters();
};

} // namespace opensn
