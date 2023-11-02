#include "LinearGraphPartitioner.h"

#include "ChiObjectFactory.h"
#include "utils/chi_utils.h"

#include "chi_log.h"

#include <cmath>

namespace chi
{

RegisterChiObject(chi, LinearGraphPartitioner);

InputParameters LinearGraphPartitioner::GetInputParameters()
{
  InputParameters params = GraphPartitioner::GetInputParameters();

  // clang-format off
  params.SetGeneralDescription("Basic linear partitioning. "
"This type of partitioner works basically only for testing. Orthogonal meshes"
" can produce decent partitioning but for unstructured grids it can be pretty"
" bad. It partitions cells based on their linear index \"global_id\" instead "
"of actually working with the graph.");
  // clang-format on
  params.SetDocGroup("Graphs");

  params.AddOptionalParameter("all_to_rank",
                              -1,
                              "If non-zero will restrict all cells to this rank, "
                              "essentially transforming this partitioner into a "
                              "single-rank partitioner.");

  return params;
}

LinearGraphPartitioner::LinearGraphPartitioner(const InputParameters& params)
  : GraphPartitioner(params),
    all_to_rank_(params.GetParamValue<int>("all_to_rank"))
{
}

/**Given a graph. Returns the partition ids of each row in the graph.*/
std::vector<int64_t> LinearGraphPartitioner::Partition(
  const std::vector<std::vector<uint64_t>>& graph,
  const std::vector<chi_mesh::Vector3>&,
  const int number_of_parts)
{
  Chi::log.Log0Verbose1() << "Partitioning with LinearGraphPartitioner";

  const std::vector<chi::SubSetInfo> sub_sets =
    chi::MakeSubSets(graph.size(), number_of_parts);

  std::vector<int64_t> pids(graph.size(), 0);

  if (all_to_rank_ < 0)
  {
    size_t n = 0;
    for (int k = 0; k < number_of_parts; ++k)
      for (size_t m = 0; m < sub_sets[k].ss_size; ++m)
        pids[n++] = k;
  }
  else
    pids.assign(graph.size(), all_to_rank_);

  Chi::log.Log0Verbose1() << "Done partitioning with LinearGraphPartitioner";
  return pids;
}

} // namespace chi