// SPDX-FileCopyrightText: 2025 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once
#include "framework/materials/multi_group_xs/multi_group_xs.h"
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
              int num_groups,
              int num_moments,
              int num_faces,
              const MultiGroupXS& xs_mapping,
              double volume,
              const std::vector<bool>& face_local_flags,
              const std::vector<int>& face_locality,
              const std::vector<const Cell*>& neighbor_cell_ptrs,
              bool cell_on_boundary)
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
    if (cell_on_boundary)
      outflow_.resize(num_faces);
  }

  size_t MapDOF(uint64_t node, uint64_t moment, uint64_t grp) const
  {
    return phi_address_ + node * num_grps_moms_ + num_groups_ * moment + grp;
  }

  const MultiGroupXS& GetXS() const { return *xs_; }

  bool IsFaceLocal(std::size_t f) const { return face_local_flags_[f]; }

  int FaceLocality(std::size_t f) const { return face_locality_[f]; }

  const Cell* FaceNeighbor(std::size_t f) const { return neighbor_cell_ptrs_[f]; }

  int GetNumNodes() const { return num_nodes_; }

  double GetVolume() const { return volume_; }

  void ZeroOutflow(std::size_t f, std::size_t g)
  {
    if (f < outflow_.size() and g < outflow_[f].size())
      outflow_[f][g] = 0.0;
  }

  void AddOutflow(std::size_t f, std::size_t g, double intS_mu_psi)
  {
    if (f < outflow_.size())
    {
      if (outflow_[f].empty())
        outflow_[f].resize(num_groups_, 0.0);
      if (g < outflow_[f].size())
        outflow_[f][g] += intS_mu_psi;
    }
  }

  double GetOutflow(std::size_t f, std::size_t g) const
  {
    if (f < outflow_.size() and g < outflow_[f].size())
      return outflow_[f][g];
    return 0.0;
  }

  void ReassignXS(const MultiGroupXS& xs) { xs_ = &xs; }

private:
  size_t phi_address_;
  int num_nodes_;
  int num_groups_;
  int num_grps_moms_;
  const MultiGroupXS* xs_;
  double volume_;
  const std::vector<bool> face_local_flags_;
  const std::vector<int> face_locality_;
  const std::vector<const Cell*> neighbor_cell_ptrs_;
  std::vector<std::vector<double>> outflow_;
};

} // namespace opensn
