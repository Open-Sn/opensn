// SPDX-FileCopyrightText: 2026 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "modules/linear_boltzmann_solvers/lbs_problem/flux_moment_layout.h"
#include "framework/utils/error.h"
#include <string>

namespace opensn
{

FluxMomentLayout::FluxMomentLayout(std::uint64_t num_local_nodes,
                                   unsigned int num_moments,
                                   unsigned int num_groups,
                                   const std::vector<GroupRange>& groupsets)
  : num_local_nodes_(num_local_nodes),
    num_moments_(num_moments),
    num_groups_(num_groups),
    num_local_dofs_(static_cast<std::size_t>(num_local_nodes) * num_moments * num_groups)
{
  OpenSnLogicalErrorIf(groupsets.empty(), "FluxMomentLayout: at least one groupset is required.");
  unsigned int next_group = 0;
  for (std::size_t gs = 0; gs < groupsets.size(); ++gs)
  {
    const auto& range = groupsets[gs];
    OpenSnLogicalErrorIf(range.first_group != next_group or range.num_groups == 0 or
                           range.num_groups > num_groups - range.first_group,
                         "FluxMomentLayout: groupset " + std::to_string(gs) + " with groups [" +
                           std::to_string(range.first_group) + ", " +
                           std::to_string(range.first_group + range.num_groups) +
                           ") does not continue a consecutive partition of " +
                           std::to_string(num_groups) + " groups.");
    next_group += range.num_groups;
  }
  OpenSnLogicalErrorIf(next_group != num_groups,
                       "FluxMomentLayout: the groupsets cover " + std::to_string(next_group) +
                         " of " + std::to_string(num_groups) + " groups.");

  // [groupset][node][moment][group in groupset]: each groupset is one block of
  // nodes x moments x (its groups) entries, and a group's stride is the size of its groupset.
  group_base_.resize(num_groups);
  group_stride_.resize(num_groups);
  group_block_.resize(num_groups);
  group_block_base_.resize(num_groups);
  block_offset_.assign(1, 0);
  for (std::size_t gs = 0; gs < groupsets.size(); ++gs)
  {
    const auto& range = groupsets[gs];
    const std::size_t block_offset = block_offset_.back();
    for (unsigned int k = 0; k < range.num_groups; ++k)
    {
      group_base_[range.first_group + k] = block_offset + k;
      group_stride_[range.first_group + k] = range.num_groups;
      group_block_[range.first_group + k] = gs;
      group_block_base_[range.first_group + k] = k;
    }
    block_offset_.push_back(block_offset + static_cast<std::size_t>(num_local_nodes) * num_moments *
                                             range.num_groups);
  }

  for (unsigned int g = 0; g < num_groups; ++g)
    node_moment_group_order_ =
      node_moment_group_order_ and group_base_[g] == g and group_stride_[g] == num_groups;
}

std::vector<std::size_t>
FluxMomentLayout::GetBlockSizes() const
{
  std::vector<std::size_t> sizes(GetNumBlocks());
  for (std::size_t b = 0; b < sizes.size(); ++b)
    sizes[b] = GetBlockSize(b);
  return sizes;
}

std::vector<double>
FluxMomentLayout::ToNodeMomentGroupOrder(std::span<const double> values) const
{
  OpenSnInvalidArgumentIf(values.size() != num_local_dofs_,
                          "FluxMomentLayout::ToNodeMomentGroupOrder: got " +
                            std::to_string(values.size()) + " values, expected " +
                            std::to_string(num_local_dofs_) + ".");
  std::vector<double> ordered(num_local_dofs_);
  std::size_t k = 0;
  for (std::uint64_t node = 0; node < num_local_nodes_; ++node)
    for (unsigned int m = 0; m < num_moments_; ++m)
      for (unsigned int g = 0; g < num_groups_; ++g)
        ordered[k++] = values[MapDOF(node, m, g)];
  return ordered;
}

std::vector<double>
FluxMomentLayout::FromNodeMomentGroupOrder(std::span<const double> values) const
{
  OpenSnInvalidArgumentIf(values.size() != num_local_dofs_,
                          "FluxMomentLayout::FromNodeMomentGroupOrder: got " +
                            std::to_string(values.size()) + " values, expected " +
                            std::to_string(num_local_dofs_) + ".");
  std::vector<double> mapped(num_local_dofs_);
  std::size_t k = 0;
  for (std::uint64_t node = 0; node < num_local_nodes_; ++node)
    for (unsigned int m = 0; m < num_moments_; ++m)
      for (unsigned int g = 0; g < num_groups_; ++g)
        mapped[MapDOF(node, m, g)] = values[k++];
  return mapped;
}

} // namespace opensn
