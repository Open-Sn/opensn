#pragma once

#include "GraphPartitioner.h"

#include "array"

namespace chi
{

class KBAGraphPartitioner : public GraphPartitioner
{
public:
  static InputParameters GetInputParameters();
  explicit KBAGraphPartitioner(const InputParameters& params);

  std::vector<int64_t> Partition(const std::vector<std::vector<uint64_t>>& graph,
                                 const std::vector<chi_mesh::Vector3>& centroids,
                                 int number_of_parts) override;

protected:
  const size_t nx_, ny_, nz_;
  const std::vector<double> xcuts_, ycuts_, zcuts_;

  struct CoordinateInfo
  {
    const std::vector<double>* cuts_;
    const size_t n_;
    const std::string coordinate_name_;
  };
  std::array<CoordinateInfo, 3> coordinate_infos_;
};

} // namespace chi


