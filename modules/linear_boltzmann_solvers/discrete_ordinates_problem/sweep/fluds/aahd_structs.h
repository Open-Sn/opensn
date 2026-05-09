// SPDX-FileCopyrightText: 2025 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep/fluds/fluds_structs.h"
#include "framework/mesh/cell/cell.h"
#include <set>
#include <tuple>

namespace opensn
{

/**
 * 64-bit integer encoding the index and the bank information of a face node.
 * \note The index must follow strict mutual-exclusion rules:
 *   - Always check boundary status BEFORE local status.
 *   - Boundary and delayed status are mutually exclusive. Both cannot be true at the same time.
 *   - Boundary and non-local status are mutually exclusive. Both cannot be true at the same time.
 *   - Non-local outgoing and delayed status are mutually exclusive. The sweep scheduler does not
 *     distinguish between delayed and non-delayed non-local outgoing faces. For non-local outgoing
 *     faces, delayed status is always set to false.
 */
class AAHD_NodeIndex : public NodeIndex
{
public:
  /// Default constructor.
  constexpr AAHD_NodeIndex() = default;

  /// Direct assign core value.
  constexpr AAHD_NodeIndex(const std::uint64_t& value) : NodeIndex(value) {}

  /**
   * Construct a non-boundary node index.
   * \param index Index into the corresponding bank. Cannot exceed 2^60 - 1.
   * \param is_outgoing Flag indicating if the node corresponds to an outgoing face.
   * \param is_local_or_reflecting Flag indicating if the index is in a local bank or
   * reflecting boundary.
   * \param is_delayed_or_angle_dependent Flag indicating if the index is in a delayed bank or an
   * angle dependent boundary.
   * \param loc_index Index into delayed/non-delayed location dependencies/successors vector for
   * non-local nodes. Cannot exceed 2^21 - 1.
   */
  AAHD_NodeIndex(std::uint64_t index,
                 bool is_outgoing,
                 bool is_boundary,
                 bool is_local_or_reflecting,
                 bool is_delayed_or_angle_dependent)
  {
    if (index >= (std::uint64_t(1) << 60) - 1)
      throw std::runtime_error("Cannot hold an index greater than 2^60.");
    SetInOut(is_outgoing);
    SetBoundary(is_boundary);
    SetLocalOrReflecting(is_local_or_reflecting);
    SetDelayedOrAngleDependent(is_delayed_or_angle_dependent);
    SetIndex(index);
  }

  /// Check if the current index corresponds to a delayed bank.
  constexpr bool IsDelayed() const noexcept
  {
    return (value_ & delayed_or_angle_dependent_bit_mask) != 0;
  }
  /// Check if the current index corresponds to an angle dependent boundary.
  constexpr bool IsAngleDependent() const noexcept { return IsDelayed(); }

  /// Check if the current index corresponds to a local bank.
  constexpr bool IsLocal() const noexcept { return (value_ & local_or_reflecting_bit_mask) != 0; }
  /// Check if the current index corresponds to a reflecting boundary.
  constexpr bool IsReflecting() const noexcept { return IsLocal(); }

  /// Get the index into the bank.
  constexpr std::uint64_t GetIndex() const noexcept { return value_ & index_bit_mask; }

private:
  /// \name Delayed/angle-dependent bit
  /// \{
  /// Third bit mask (``001`` followed by 61 zeros) - Bit 61
  static constexpr std::uint64_t delayed_or_angle_dependent_bit_mask = std::uint64_t(1) << (64 - 3);
  /// Encode the value as delayed or angle dependent boundary.
  constexpr void SetDelayedOrAngleDependent(bool is_delayed_or_angle_dependent) noexcept
  {
    if (is_delayed_or_angle_dependent)
      value_ |= delayed_or_angle_dependent_bit_mask;
    else
      value_ &= ~delayed_or_angle_dependent_bit_mask;
  }
  /// \}

  /// \name Local/reflecting bit
  /// \{
  /// Fourth bit mask (``0001`` followed by 60 zeros) - Bit 60
  static constexpr std::uint64_t local_or_reflecting_bit_mask = std::uint64_t(1) << (64 - 4);
  /// Encode the value as local.
  constexpr void SetLocalOrReflecting(bool is_local_or_reflecting) noexcept
  {
    if (is_local_or_reflecting)
      value_ |= local_or_reflecting_bit_mask;
    else
      value_ &= ~local_or_reflecting_bit_mask;
  }
  /// \}

  /// \name Index bits
  /// \{
  /// Index bit mask (``1`` at the last 60 bits).
  static constexpr std::uint64_t index_bit_mask = (std::uint64_t(1) << (64 - 4)) - 1;
  /// Encode the index.
  constexpr void SetIndex(const std::uint64_t& index) noexcept
  {
    value_ &= ~index_bit_mask;
    value_ |= (index & index_bit_mask);
  }
  /// \}
};

/// Node on a face (directed edge) connecting an upwind and downwind cell.
struct AAHD_DirectedEdgeNode
{
  /// Upwind face node.
  FaceNode upwind_node;
  /// Downwind face node.
  FaceNode downwind_node;

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
  FaceNode node;

private:
  /// Tuple for ordering the non-local face nodes.
  std::tuple<std::uint64_t, std::uint64_t, std::uint64_t> global_ordering_;
};

/**
 * Set of pointers for AAHD FLUDS.
 */
struct AAHD_FLUDSPointerSet : public FLUDSPointerSet
{
  /// Pointer to delayed local angular fluxes.
  double* __restrict__ delayed_local_psi = nullptr;
  /// Pointer to delayed old local angular fluxes.
  double* __restrict__ delayed_local_psi_old = nullptr;
  /// Pointer to non-local old delayed incoming angular fluxes.
  double* __restrict__ nonlocal_delayed_incoming_psi_old = nullptr;

  /// Get pointer to the incoming angular flux (if the face is not incoming, a nullptr is returned).
  constexpr double* GetIncomingFluxPointer(const AAHD_NodeIndex& node_index,
                                           unsigned int angle_group_idx,
                                           unsigned int group_idx,
                                           double* __restrict__ boundary,
                                           const std::uint64_t* __restrict__ boundary_offset,
                                           bool is_surface_source_active) const noexcept
  {
    // boundary case
    if (node_index.IsBoundary())
    {
      if (node_index.IsReflecting() or is_surface_source_active)
      {
        unsigned int location = (node_index.IsAngleDependent()) ? angle_group_idx : group_idx;
        return boundary + boundary_offset[node_index.GetIndex()] + location;
      }
      return boundary + group_idx;
    }
    // local case
    if (node_index.IsLocal())
    {
      if (node_index.IsDelayed())
      {
        return delayed_local_psi_old + node_index.GetIndex() * stride_size + angle_group_idx;
      }
      else
      {
        return local_psi + node_index.GetIndex() * stride_size + angle_group_idx;
      }
    }
    // non-local case
    else
    {
      if (node_index.IsDelayed())
      {
        return nonlocal_delayed_incoming_psi_old + node_index.GetIndex() * stride_size +
               angle_group_idx;
      }
      else
      {
        return nonlocal_incoming_psi + node_index.GetIndex() * stride_size + angle_group_idx;
      }
    }
  }

  /// Get pointer to the outgoing angular flux (if the face is not outgoing, a nullptr is returned).
  constexpr double*
  GetOutgoingFluxPointer(const AAHD_NodeIndex& node_index,
                         unsigned int angle_group_idx,
                         double* __restrict__ boundary,
                         const std::uint64_t* __restrict__ boundary_offset) const noexcept
  {
    // boundary case
    if (node_index.IsBoundary())
    {
      return boundary + boundary_offset[node_index.GetIndex()] + angle_group_idx;
    }
    // local case
    if (node_index.IsLocal())
    {
      if (node_index.IsDelayed())
      {
        return delayed_local_psi + node_index.GetIndex() * stride_size + angle_group_idx;
      }
      else
      {
        return local_psi + node_index.GetIndex() * stride_size + angle_group_idx;
      }
    }
    // non-local case
    else
    {
      return nonlocal_outgoing_psi + node_index.GetIndex() * stride_size + angle_group_idx;
    }
  }
};

} // namespace opensn

namespace std
{

/// Print non-local face node.
ostream& operator<<(ostream& out, const opensn::AAHD_NonLocalFaceNode& n);

/// Print node index.
ostream& operator<<(ostream& out, const opensn::AAHD_NodeIndex& e);

} // namespace std
