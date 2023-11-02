#include "chi_utils.h"

#include <cmath>

#define scdouble static_cast<double>

namespace chi
{

std::vector<SubSetInfo>
MakeSubSets(size_t num_items, size_t desired_num_subsets)
{
  std::vector<SubSetInfo> ss_infos;
  const std::size_t div = std::floor(scdouble(num_items) / scdouble(desired_num_subsets));
  const std::size_t rem = num_items % desired_num_subsets;

  for (size_t i = 0; i < desired_num_subsets; ++i)
    ss_infos.push_back({0, 0, div});
  for (size_t j = 0; j < rem; ++j)
    ss_infos[j].ss_size += 1;

  size_t check_sum = 0;
  for (size_t i = 0; i < desired_num_subsets; ++i)
  {
    ss_infos[i].ss_begin = check_sum;
    check_sum += ss_infos[i].ss_size;
    ss_infos[i].ss_end = check_sum - 1;
  }

  return ss_infos;
}

} // namespace chi