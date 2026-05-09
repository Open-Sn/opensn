// SPDX-FileCopyrightText: 2025 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once
#include "framework/materials/multi_group_xs/multi_group_xs.h"
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
  CellLBSView(size_t phi_address,
              int num_nodes,
              unsigned int num_groups,
              unsigned int num_moments,
              const MultiGroupXS& xs_mapping,
              double volume,
              const std::vector<bool>& face_local_flags,
              const std::vector<int>& face_locality,
              const std::vector<const Cell*>& neighbor_cell_ptrs)
    : phi_address_(phi_address),
      num_nodes_(num_nodes),
      num_groups_(num_groups),
      num_grps_moms_(num_groups * num_moments),
      xs_(&xs_mapping),
      volume_(volume),
      face_local_flags_(face_local_flags),
      face_locality_(face_locality),
      neighbor_cell_ptrs_(neighbor_cell_ptrs)
  {
  }

  size_t MapDOF(uint64_t node, uint64_t moment, unsigned int grp) const
  {
    return phi_address_ + node * num_grps_moms_ + num_groups_ * moment + grp;
  }

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
  size_t phi_address_;
  int num_nodes_;
  unsigned int num_groups_;
  unsigned int num_grps_moms_;
  const MultiGroupXS* xs_;
  double volume_;
  const std::vector<bool> face_local_flags_;
  const std::vector<int> face_locality_;
  const std::vector<const Cell*> neighbor_cell_ptrs_;
};

} // namespace opensn
