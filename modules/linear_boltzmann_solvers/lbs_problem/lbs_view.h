// SPDX-FileCopyrightText: 2025 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once
#include "framework/materials/multi_group_xs/multi_group_xs.h"
#include "modules/linear_boltzmann_solvers/lbs_problem/flux_moment_layout.h"
#include <cassert>
#include <cstddef>
#include <cstdint>
#include <vector>

namespace opensn
{

class Cell;

/// Transport view of a cell.
class CellLBSView
{
public:
  /**
   * `layout` is the problem's flux-moment map and must outlive the view. `first_node` is the local
   * node index of the cell's node 0.
   */
  CellLBSView(const FluxMomentLayout& layout,
              std::uint64_t first_node,
              int num_nodes,
              const MultiGroupXS& xs_mapping,
              double volume,
              const std::vector<bool>& face_local_flags,
              const std::vector<int>& face_locality,
              const std::vector<const Cell*>& neighbor_cell_ptrs)
    : layout_(&layout),
      first_node_(first_node),
      num_nodes_(num_nodes),
      xs_(&xs_mapping),
      volume_(volume),
      face_local_flags_(face_local_flags),
      face_locality_(face_locality),
      neighbor_cell_ptrs_(neighbor_cell_ptrs)
  {
  }

  /// Index into the flux-moment vectors of cell node `node`, moment `moment`, and group `grp`.
  size_t MapDOF(uint64_t node, uint64_t moment, uint64_t grp) const
  {
    assert(node < static_cast<uint64_t>(num_nodes_) && "CellLBSView::MapDOF node out of range.");
    return layout_->MapDOF(first_node_ + node, moment, grp);
  }

  /// Local node index of the cell's node 0 in the flux-moment layout.
  std::uint64_t GetFirstNode() const { return first_node_; }

  const MultiGroupXS& GetXS() const { return *xs_; }

  bool IsFaceLocal(std::size_t f) const
  {
    assert(f < face_local_flags_.size() && "CellLBSView::IsFaceLocal face index out of range.");
    return face_local_flags_[f];
  }

  int FaceLocality(std::size_t f) const
  {
    assert(f < face_locality_.size() && "CellLBSView::FaceLocality face index out of range.");
    return face_locality_[f];
  }

  const Cell* FaceNeighbor(std::size_t f) const
  {
    assert(f < neighbor_cell_ptrs_.size() && "CellLBSView::FaceNeighbor face index out of range.");
    return neighbor_cell_ptrs_[f];
  }

  int GetNumNodes() const { return num_nodes_; }

  double GetVolume() const { return volume_; }

  void ReassignXS(const MultiGroupXS& xs) { xs_ = &xs; }

private:
  const FluxMomentLayout* layout_;
  std::uint64_t first_node_;
  int num_nodes_;
  const MultiGroupXS* xs_;
  double volume_;
  const std::vector<bool> face_local_flags_;
  const std::vector<int> face_locality_;
  const std::vector<const Cell*> neighbor_cell_ptrs_;
};

} // namespace opensn
