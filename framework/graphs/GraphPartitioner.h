#pragma once

#include "ChiObject.h"

namespace chi_mesh
{
struct Vector3;
}

namespace chi
{

/**Abstract base class for all partitioners*/
class GraphPartitioner : public ChiObject
{
public:
  /**Given a graph. Returns the partition ids of each row in the graph.*/
  virtual std::vector<int64_t> Partition(const std::vector<std::vector<uint64_t>>& graph,
                                         const std::vector<chi_mesh::Vector3>& centroids,
                                         int number_of_parts) = 0;

protected:
  static InputParameters GetInputParameters();
  explicit GraphPartitioner(const InputParameters& params);
};

} // namespace chi
