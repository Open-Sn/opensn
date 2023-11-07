#pragma once

#include "framework/graphs/graph_partitioner.h"

namespace chi
{

class LinearGraphPartitioner : public GraphPartitioner
{
public:
  static InputParameters GetInputParameters();
  explicit LinearGraphPartitioner(const InputParameters& params);

  std::vector<int64_t> Partition(const std::vector<std::vector<uint64_t>>& graph,
                                 const std::vector<chi_mesh::Vector3>& centroids,
                                 int number_of_parts) override;

protected:
  const int all_to_rank_;
};

} // namespace chi
