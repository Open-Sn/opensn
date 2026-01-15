// SPDX-FileCopyrightText: 2025 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep/fluds/device_structs.h"
#include <set>
#include <tuple>

namespace opensn
{

using AAHD_FaceNode = DeviceFaceNode;

using AAHD_NodeIndex = DeviceNodeIndex;

/// Node on a face (directed edge) connecting an upwind and downwind cell.
struct AAHD_DirectedEdgeNode
{
  /// Upwind face node.
  AAHD_FaceNode upwind_node;
  /// Downwind face node.
  AAHD_FaceNode downwind_node;

  /// Check if the directed edge is initialized.
  constexpr bool IsInitialized() const
  {
    return upwind_node.IsInitialized() && downwind_node.IsInitialized();
  }
};

/// Class representing non-local face node between 2 cells in different partitions.
class AAHD_NonLocalFaceNode
{
public:
  /// Default constructor.
  AAHD_NonLocalFaceNode() = default;

  /// Member constructor.
  AAHD_NonLocalFaceNode(std::uint64_t cell_global_idx,
                        std::uint32_t cell_local_idx,
                        std::uint16_t face_idx,
                        std::uint16_t face_node_idx,
                        std::uint64_t neighbor_global_idx,
                        std::uint16_t neighbor_face_idx,
                        std::uint16_t neighbor_face_node_idx)
    : node(cell_local_idx, face_idx, face_node_idx)
  {
    std::uint32_t packed_cell_idx = (std::uint32_t(face_idx) << 16) | face_node_idx;
    std::uint32_t packed_neighbor_idx =
      (std::uint32_t(neighbor_face_idx) << 16) | neighbor_face_node_idx;
    if (cell_global_idx < neighbor_global_idx)
    {
      global_ordering_ =
        std::make_tuple(cell_global_idx,
                        neighbor_global_idx,
                        (std::uint64_t(packed_cell_idx) << 32) | packed_neighbor_idx);
    }
    else if (neighbor_global_idx < cell_global_idx)
    {
      global_ordering_ =
        std::make_tuple(neighbor_global_idx,
                        cell_global_idx,
                        (std::uint64_t(packed_neighbor_idx) << 32) | packed_cell_idx);
    }
    else
      throw std::runtime_error("Non local face node requires cell with different global index.\n");
  }

  /// Comparison operator for ordering.
  bool operator<(const AAHD_NonLocalFaceNode& other) const
  {
    return global_ordering_ < other.global_ordering_;
  }

  /// Get the global ordering tuple.
  const std::tuple<std::uint64_t, std::uint64_t, std::uint64_t>& GetGlobalOrdering() const
  {
    return global_ordering_;
  }

  /// Corresponding face node.
  AAHD_FaceNode node;

private:
  /// Tuple for ordering the non-local face nodes.
  std::tuple<std::uint64_t, std::uint64_t, std::uint64_t> global_ordering_;
};

/**
 * \brief Set of pointers for AAHD FLUDS.
 */
struct AAHD_FLUDSPointerSet
{
  /// Pointer to local angular fluxes.
  double* local_psi = nullptr;
  /// Pointer to delayed local angular fluxes.
  double* delayed_local_psi = nullptr;
  /// Pointer to delayed old local angular fluxes.
  double* delayed_local_psi_old = nullptr;
  /// Pointer to non-local incoming angular fluxes.
  double* nonlocal_incoming_psi = nullptr;
  /// Pointer to non-local delayed incoming angular fluxes.
  double* nonlocal_delayed_incoming_psi = nullptr;
  /// Pointer to non-local old delayed incoming angular fluxes.
  double* nonlocal_delayed_incoming_psi_old = nullptr;
  /// Pointer to non-local outgoing angular fluxes.
  double* nonlocal_outgoing_psi = nullptr;
  /// Pointer to boundary angular fluxes.
  double* boundary_psi = nullptr;
  /// Stride size (number of groups times number of angles).
  std::uint64_t stride_size = 0;

  /// Get pointer to the incoming angular flux (if the face is not incoming, a nullptr is returned).
  constexpr double* GetIncomingFluxPointer(const AAHD_NodeIndex& node_index) const noexcept
  {
    // undefined case (parallel face) : all tests are true
    if (node_index.IsUndefined())
    {
      return nullptr;
    }
    // outgoing case : nullptr
    if (node_index.IsOutgoing())
    {
      return nullptr;
    }
    // boundary case : non-delayed and local tests are true
    if (node_index.IsBoundary())
    {
      return boundary_psi + node_index.GetIndex() * stride_size;
    }
    // local case
    if (node_index.IsLocal())
    {
      if (node_index.IsDelayed())
      {
        return delayed_local_psi_old + node_index.GetIndex() * stride_size;
      }
      else
      {
        return local_psi + node_index.GetIndex() * stride_size;
      }
    }
    // non-local case
    else
    {
      if (node_index.IsDelayed())
      {
        return nonlocal_delayed_incoming_psi_old + node_index.GetIndex() * stride_size;
      }
      else
      {
        return nonlocal_incoming_psi + node_index.GetIndex() * stride_size;
      }
    }
  }

  /// Get pointer to the outgoing angular flux (if the face is not incoming, a nullptr is returned).
  constexpr double* GetOutgoingFluxPointer(const AAHD_NodeIndex& node_index) const noexcept
  {
    // undefined case (parallel face) : all tests are true
    if (node_index.IsUndefined())
    {
      return nullptr;
    }
    // incoming case : nullptr
    if (!node_index.IsOutgoing())
    {
      return nullptr;
    }
    // boundary case : non-delayed and local tests are true
    if (node_index.IsBoundary())
    {
      return boundary_psi + node_index.GetIndex() * stride_size;
    }
    // local case
    if (node_index.IsLocal())
    {
      if (node_index.IsDelayed())
      {
        return delayed_local_psi + node_index.GetIndex() * stride_size;
      }
      else
      {
        return local_psi + node_index.GetIndex() * stride_size;
      }
    }
    // non-local case
    else
    {
      return nonlocal_outgoing_psi + node_index.GetIndex() * stride_size;
    }
  }
};

} // namespace opensn

namespace std
{

/// Print non-local face node.
ostream& operator<<(ostream& out, const opensn::AAHD_NonLocalFaceNode& n);

} // namespace std
