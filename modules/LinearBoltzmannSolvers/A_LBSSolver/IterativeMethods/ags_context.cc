#include "ags_context.h"

#include "A_LBSSolver/lbs_solver.h"
#include "wgs_context.h"

#include <petscksp.h>

#define GetGSContextPtr(x) std::dynamic_pointer_cast<WGSContext<Mat, Vec, KSP>>(x)

namespace lbs
{

template <>
std::pair<int64_t, int64_t>
AGSContext<Mat, Vec, KSP>::SystemSize()
{
  const size_t local_node_count = lbs_solver_.LocalNodeCount();
  const size_t globl_node_count = lbs_solver_.GlobalNodeCount();

  std::set<int> groupset_list_group_ids;
  for (auto& wgs_solver : sub_solvers_list_)
  {
    auto gs_context_ptr = GetGSContextPtr(wgs_solver->GetContext());
    for (const auto& group : gs_context_ptr->groupset_.groups_)
      groupset_list_group_ids.insert(group.id_);
  }

  const size_t gslist_num_groups = groupset_list_group_ids.size();

  const size_t num_moments = lbs_solver_.NumMoments();

  const size_t local_size = local_node_count * gslist_num_groups * num_moments;
  const size_t globl_size = globl_node_count * gslist_num_groups * num_moments;

  return {static_cast<int64_t>(local_size), static_cast<int64_t>(globl_size)};
}

template <>
void
AGSContext<Mat, Vec, KSP>::SetPreconditioner(KSP& solver)
{
}

template <>
int
AGSContext<Mat, Vec, KSP>::MatrixAction(Mat& matrix, Vec& vector, Vec& action)
{

  return 0;
}

} // namespace lbs
