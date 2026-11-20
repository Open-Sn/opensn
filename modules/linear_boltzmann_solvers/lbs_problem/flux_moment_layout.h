// SPDX-FileCopyrightText: 2026 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include <cassert>
#include <cstddef>
#include <cstdint>
#include <vector>

namespace opensn
{

/**
 * Index map of the local flux-moment vectors: the scalar flux moments, the source moments, the
 * external source moments, the uncollided flux moments, and any temporary with the same size.
 *
 * The map takes a local node index, a moment index, and a group index. It is the only description
 * of the storage order, so code that indexes these vectors must go through it rather than assume
 * an order.
 *
 * The map guarantees one property: at a fixed node and moment, consecutive groups of one groupset
 * are stored at consecutive indices. A pointer to the first group of a groupset therefore addresses
 * every group of that groupset. Groups of different groupsets have no defined relative position.
 *
 * The current storage order is [node][moment][group].
 */
class FluxMomentLayout
{
public:
  /// The groups of one groupset.
  struct GroupRange
  {
    unsigned int first_group = 0;
    unsigned int num_groups = 0;
  };

  FluxMomentLayout() = default;

  /**
   * Builds the map for `num_local_nodes` nodes, `num_moments` moments, and `num_groups` groups.
   * `groupsets` must partition the groups consecutively starting at group 0, in order; otherwise
   * throws std::logic_error.
   */
  FluxMomentLayout(std::uint64_t num_local_nodes,
                   unsigned int num_moments,
                   unsigned int num_groups,
                   const std::vector<GroupRange>& groupsets);

  /// Number of local entries: nodes x moments x groups.
  std::size_t GetNumLocalDOFs() const { return num_local_dofs_; }

  std::uint64_t GetNumLocalNodes() const { return num_local_nodes_; }

  unsigned int GetNumMoments() const { return num_moments_; }

  unsigned int GetNumGroups() const { return num_groups_; }

  /// Index of local node `node`, moment `moment`, and group `group`.
  std::size_t MapDOF(std::uint64_t node, std::uint64_t moment, std::uint64_t group) const
  {
    assert(node < num_local_nodes_ && "FluxMomentLayout::MapDOF node index out of range.");
    assert(moment < num_moments_ && "FluxMomentLayout::MapDOF moment index out of range.");
    assert(group < num_groups_ && "FluxMomentLayout::MapDOF group index out of range.");
    return group_base_[group] + (node * num_moments_ + moment) * group_stride_[group];
  }

  /**
   * Distance between the indices of consecutive moments of `group` at one node. The distance
   * between consecutive nodes is this value times the number of moments.
   */
  std::size_t GetMomentStride(std::uint64_t group) const
  {
    assert(group < num_groups_ && "FluxMomentLayout::GetMomentStride group index out of range.");
    return group_stride_[group];
  }

  /**
   * True when the layout is [node][moment][group] order, the order of the flux-moment unknown
   * manager. Code that indexes these vectors through that unknown manager, such as PWLD ghost
   * communication, requires it.
   */
  bool IsNodeMomentGroupOrder() const { return node_moment_group_order_; }

  /**
   * Returns a copy of `values`, given in this layout, in [node][moment][group] order: the order of
   * flux-moment files and of the Python flux accessors. Throws std::invalid_argument when
   * `values` does not have GetNumLocalDOFs() entries.
   */
  std::vector<double> ToNodeMomentGroupOrder(const std::vector<double>& values) const;

  /**
   * Returns a copy of `values`, given in [node][moment][group] order, in this layout. Throws
   * std::invalid_argument when `values` does not have GetNumLocalDOFs() entries.
   */
  std::vector<double> FromNodeMomentGroupOrder(const std::vector<double>& values) const;

private:
  std::uint64_t num_local_nodes_ = 0;
  unsigned int num_moments_ = 0;
  unsigned int num_groups_ = 0;
  std::size_t num_local_dofs_ = 0;
  bool node_moment_group_order_ = true;
  /// Index of node 0, moment 0 for each group.
  std::vector<std::size_t> group_base_;
  /// Moment stride for each group.
  std::vector<std::size_t> group_stride_;
};

} // namespace opensn
