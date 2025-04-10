// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include "framework/graphs/graph_partitioner.h"
#include "framework/object_factory.h"
#include <array>
#include <memory>

namespace opensn
{

class KBAGraphPartitioner : public GraphPartitioner
{
public:
  explicit KBAGraphPartitioner(const InputParameters& params);

  std::vector<int64_t> Partition(const std::vector<std::vector<uint64_t>>& graph,
                                 const std::vector<Vector3>& centroids,
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

public:
  static InputParameters GetInputParameters();
  static std::shared_ptr<KBAGraphPartitioner> Create(const ParameterBlock& params);
};

} // namespace opensn
