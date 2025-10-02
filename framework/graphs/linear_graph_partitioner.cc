// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "framework/graphs/linear_graph_partitioner.h"
#include "framework/utils/utils.h"
#include "framework/object_factory.h"
#include "framework/runtime.h"
#include "framework/logging/log.h"
#include <cmath>

namespace opensn
{

OpenSnRegisterObjectInNamespace(mesh, LinearGraphPartitioner);

InputParameters
LinearGraphPartitioner::GetInputParameters()
{
  InputParameters params = GraphPartitioner::GetInputParameters();

  params.SetGeneralDescription(
    "Basic linear partitioning. This type of partitioner works basically only for testing. "
    "Orthogonal meshes can produce decent partitioning but for unstructured grids it can be pretty "
    "bad. It partitions cells based on their linear index \"global_id\" instead of actually "
    "working with the graph.");

  params.AddOptionalParameter("all_to_rank",
                              -1,
                              "If non-zero will restrict all cells to this rank, "
                              "essentially transforming this partitioner into a "
                              "single-rank partitioner.");

  return params;
}

std::shared_ptr<LinearGraphPartitioner>
LinearGraphPartitioner::Create(const ParameterBlock& params)
{
  auto& factory = opensn::ObjectFactory::GetInstance();
  return factory.Create<LinearGraphPartitioner>("mesh::LinearGraphPartitioner", params);
}

LinearGraphPartitioner::LinearGraphPartitioner(const InputParameters& params)
  : GraphPartitioner(params), all_to_rank_(params.GetParamValue<int>("all_to_rank"))
{
}

std::vector<int>
LinearGraphPartitioner::Partition(const std::vector<std::vector<uint64_t>>& graph,
                                  const std::vector<Vector3>& /* centroids */,
                                  const int number_of_parts)
{
  log.Log0Verbose1() << "Partitioning with LinearGraphPartitioner";

  const std::vector<SubSetInfo> sub_sets = MakeSubSets(graph.size(), number_of_parts);

  std::vector<int> pids(graph.size(), 0);

  if (all_to_rank_ < 0)
  {
    size_t n = 0;
    for (int k = 0; k < number_of_parts; ++k)
      for (size_t m = 0; m < sub_sets[k].ss_size; ++m)
        pids[n++] = k;
  }
  else
    pids.assign(graph.size(), all_to_rank_);

  log.Log0Verbose1() << "Done partitioning with LinearGraphPartitioner";
  return pids;
}

} // namespace opensn
