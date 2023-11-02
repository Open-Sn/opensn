#ifndef CHITECH_PETSCGRAPHPARTITIONER_H
#define CHITECH_PETSCGRAPHPARTITIONER_H

#include "GraphPartitioner.h"

namespace chi
{

class PETScGraphPartitioner : public GraphPartitioner
{
public:
  static InputParameters GetInputParameters();
  explicit PETScGraphPartitioner(const InputParameters& params);

  std::vector<int64_t> Partition(const std::vector<std::vector<uint64_t>>& graph,
                                 const std::vector<chi_mesh::Vector3>& centroids,
                                 int number_of_parts) override;

protected:
  const std::string type_;
};

} // namespace chi

#endif // CHITECH_PETSCGRAPHPARTITIONER_H
