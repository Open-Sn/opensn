// SPDX-FileCopyrightText: 2025 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include <cstdint>
#include <limits>
#include <ostream>
#include <set>
#include <stdexcept>
#include <tuple>

namespace opensn
{

/**
 * @brief Node on a face, which is also an edge in the sweep graph.
 * @details This class representing a node on a face. It stores the cell local index, the face index
 * within the cell and the index of the node relative to the face as compact 64-bit integer. This
 * class abstract the node for usage with std::map and std::set.
 */
class AAHD_FaceNode
{
public:
  /// Default constructor.
  constexpr AAHD_FaceNode() = default;
  /// Member constructor.
  constexpr AAHD_FaceNode(const std::uint32_t& cell_idx,
                     const std::uint16_t& face_idx,
                     const std::uint16_t& face_node_idx)
  {
    value_ = 0;
    value_ |= static_cast<std::uint64_t>(cell_idx) << 32;
    value_ |= static_cast<std::uint64_t>(face_idx) << 16;
    value_ |= static_cast<std::uint64_t>(face_node_idx);
  }

  /// Comparison operator for ordering.
  constexpr bool operator<(const AAHD_FaceNode& other) const { return value_ < other.value_; }

  /// Equality operator.
  constexpr bool operator==(const AAHD_FaceNode& other) const { return value_ == other.value_; }

  /// Get cell local index.
  constexpr std::uint32_t GetCellIndex() const { return static_cast<std::uint32_t>(value_ >> 32); }

  /// Get face index.
  constexpr std::uint16_t GetFaceIndex() const
  {
    return static_cast<std::uint16_t>((value_ >> 16) & 0xFFFFu);
  }

  /// Get face node index.
  constexpr std::uint16_t GetFaceNodeIndex() const
  {
    return static_cast<std::uint16_t>(value_ & 0xFFFFu);
  }

  /// Check if the face node is initialized.
  constexpr bool IsInitialized() const
  {
    return value_ != std::numeric_limits<std::uint64_t>::max();
  }

private:
  /// Core value.
  std::uint64_t value_ = std::numeric_limits<std::uint64_t>::max();
};

/// Node on a face (directed edge) connecting an upwind and downwind cell.
struct AAHD_DirectedEdgeNode
{
  /// Upwind face node.
  AAHD_FaceNode upwind_node;
  /// Downwind face node.
  AAHD_FaceNode downwind_node;

  /// Check if the directed edge is coming from a given set of cell local indices.
  inline bool IsComingFrom(const std::set<std::uint32_t>& sorted_level) const
  {
    return sorted_level.contains(upwind_node.GetCellIndex());
  }

  /// Check if the directed edge is initialized.
  constexpr bool IsInitialized() const
  {
    return upwind_node.IsInitialized() && downwind_node.IsInitialized();
  }
};

/// Class reprent non-local face node between 2 cells in different partitions.
class AAHD_NonLocalFaceNode
{
public:
  /// Default constructor.
  AAHD_NonLocalFaceNode() = default;

  /// Member constructor.
  AAHD_NonLocalFaceNode(const std::uint64_t& cell_global_idx,
                   const std::uint32_t& cell_local_idx,
                   const std::uint16_t& face_idx,
                   const std::uint16_t& face_node_idx,
                   const std::uint64_t& neighbor_global_idx,
                   const std::uint16_t& neighbor_face_idx,
                   const std::uint16_t& neighbor_face_node_idx)
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
  inline bool operator<(const AAHD_NonLocalFaceNode& other) const
  {
    return global_ordering_ < other.global_ordering_;
  }

  /// Get the global ordering tuple.
  inline const std::tuple<std::uint64_t, std::uint64_t, std::uint64_t>& GetGlobalOrdering() const
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
 * @brief 64-bit integer encoding the index and the bank information of a face node.
 * @note The index must follow strict mutual-exclusion rules:
 *   - Always check boundary status BEFORE local status. Boundary instances are assigned
 *     as local, but their locations are actually set to an undefined value.
 *   - Boundary and delayed status are mutually exclusive. Both cannot be true at the same time.
 *   - Boundary and non-local status are mutually exclusive. Both cannot be true at the same time.
 *   - Non-local outgoing and delayed status are mutually exclusive. The sweep scheduler does not
 *     distinguish between delayed and non-delayed non-local outgoing faces. For non-local outgoing
 *     faces, delayed status is always set to false.
 */
class AAHD_NodeIndex
{
public:
  /// Integer representation of the current location (``1`` on the last 21 bits).
  static inline constexpr std::uint32_t current_loc = (std::uint32_t(1) << 21) - 1;

  /// Default constructor
  AAHD_NodeIndex() = default;

  /// Direct assign core value
  constexpr AAHD_NodeIndex(const std::uint64_t& value) : value_(value) {}

  /**
   * Create an index for non-boundary node.
   * @param index Index into the corresponding bank. Cannot exceed 2^40 - 1.
   * @param is_outgoing Flag indicating if the node corresponds to an outgoing face.
   * @param is_delayed Flag indicating if the index is in a delayed bank.
   * @param loc_index Index into delayed/non-delayed location dependencies/successors vector for
   * non-local nodes. Cannot exceed 2^21 - 1.
   */
  static inline AAHD_NodeIndex CreateIndex(const std::uint64_t& index,
                                      bool is_outgoing,
                                      bool is_delayed = false,
                                      const std::uint32_t& loc_index = current_loc)
  {
    AAHD_NodeIndex result;
    if (loc_index > current_loc)
      throw std::runtime_error("Cannot accomodate MPI with rank greater than 2^21.");
    if (index >= (std::uint64_t(1) << 40) - 1)
      throw std::runtime_error("Cannot hold an index greater than 2^40.");
    result.SetInOut(is_outgoing);
    result.SetDelayed(is_delayed);
    result.SetBoundary(false);
    result.SetLocation(loc_index);
    result.SetIndex(index);
    if (is_outgoing && loc_index != current_loc)
      result.SetDelayed(false);
    return result;
  }

  /**
   * Create an index for boundary node.
   * @param index Index into the corresponding bank. Cannot exceed 2^40 - 1.
   * @param is_outgoing Flag indicating if the node corresponds to an outgoing face.
   */
  static inline AAHD_NodeIndex CreateBoundaryIndex(const std::uint64_t& index, bool is_outgoing)
  {
    AAHD_NodeIndex result;
    result.SetInOut(is_outgoing);
    result.SetDelayed(false);
    result.SetBoundary(true);
    result.SetIndex(index);
    return result;
  }

  /// Check if the current index is in undefined form.
  constexpr bool IsUndefined() const noexcept
  {
    return value_ == std::numeric_limits<std::uint64_t>::max();
  }

  /// Check if the current index is incoming or outgoing.
  constexpr bool IsOutGoing() const noexcept { return (value_ & inout_bit_mask) != 0; }

  /// Check if the current index is corresponding to a delayed bank.
  constexpr bool IsDelayed() const noexcept { return (value_ & delayed_bit_mask) != 0; }

  /// Check if the current index is corresponding to boundary.
  constexpr bool IsBoundary() const noexcept { return (value_ & boundary_bit_mask) != 0; }

  /// Check if the current index is corresponding to a local bank.
  constexpr bool IsLocal() const noexcept { return GetRank() == current_loc; }

  /// Get the rank encoded by the index (not working for current rank).
  constexpr std::uint32_t GetRank() const noexcept { return (value_ & loc_bit_mask) >> 40; }

  /// Get the index into the bank.
  constexpr std::uint64_t GetIndex() const noexcept { return value_ & index_bit_mask; }

  /// Get the core value.
  constexpr std::uint64_t GetCoreValue() const noexcept { return value_; }

private:
  /// @name Incoming/outgoing bit
  /// @{
  /// First bit mask (``1`` followed by 63 zeros).
  static constexpr std::uint64_t inout_bit_mask = std::uint64_t(1) << (64 - 1);
  /// Encode the value as incoming or outgoing.
  constexpr void SetInOut(bool is_outgoing) noexcept
  {
    if (is_outgoing)
      value_ |= inout_bit_mask;
    else
      value_ &= ~inout_bit_mask;
  }
  /// @}

  /// @name Boundary bit
  /// @{
  /// Second bit mask (``01`` followed by 62 zeros).
  static constexpr std::uint64_t boundary_bit_mask = std::uint64_t(1) << (64 - 2);
  /// Encode the value as boundary.
  constexpr void SetBoundary(bool is_boundary) noexcept
  {
    if (is_boundary)
      value_ |= boundary_bit_mask;
    else
      value_ &= ~boundary_bit_mask;
  }
  /// @}

  /// @name Delayed bit
  /// @{
  /// Third bit mask (``001`` followed by 61 zeros).
  static constexpr std::uint64_t delayed_bit_mask = std::uint64_t(1) << (64 - 3);
  /// Encode the value as delayed.
  constexpr void SetDelayed(bool is_delayed) noexcept
  {
    if (is_delayed)
      value_ |= delayed_bit_mask;
    else
      value_ &= ~delayed_bit_mask;
  }
  /// @}

  /// @name Location bits
  /// @{
  /// Location bit mask (``1`` from the 4th to 24th bit).
  static constexpr std::uint64_t loc_bit_mask = std::uint64_t(current_loc) << (64 - 3 - 21);
  /// Encode the location.
  constexpr void SetLocation(const std::uint32_t& loc_index) noexcept
  {
    std::uint64_t v = std::uint64_t(loc_index & current_loc) << (64 - 3 - 21);
    value_ &= ~loc_bit_mask;
    value_ |= v;
  }
  /// @}

  /// @name Index bits
  /// @{
  /// Index bit mask (``1`` at the last 40 bits)
  static constexpr std::uint64_t index_bit_mask = (std::uint64_t(1) << (64 - 3 - 21)) - 1;
  /// Encode the index.
  constexpr void SetIndex(const std::uint64_t& index) noexcept
  {
    value_ &= ~index_bit_mask;
    value_ |= (index & index_bit_mask);
  }
  /// @}

  /**
   * @brief Core value.
   * @details 1 bit for incoming/outgoing; 1 bit for delayed or non delayed; 1 bit for boundary or
   * not; 21 bit (~2.1M) for MPI rank (all ones if local); 40 bits (~4T) for index.
   */
  std::uint64_t value_ = std::numeric_limits<std::uint64_t>::max();
};

/**
 * @brief Set of pointers for AAHD FLUDS.
 */
struct AAHD_FLUDSPointerSet
{
  /// Pointer to local angular fluxes.
  double* local_psi = nullptr;
  /// Pointer to delayed local angular fluxes.
  double* delayed_local_psi = nullptr;
  /// Pointer to delayed old local angular fluxes.
  double* delayed_local_psi_old = nullptr;
  /// Pointer to non-local incoming offset.
  std::uint64_t* nonlocal_incoming_offsets = nullptr;
  /// Pointer to non-local incoming angular fluxes.
  double* nonlocal_incoming_psi = nullptr;
  /// Pointer to non-local delayed incoming offset.
  std::uint64_t* nonlocal_delayed_incoming_offsets = nullptr;
  /// Pointer to non-local delayed incoming angular fluxes.
  double* nonlocal_delayed_incoming_psi = nullptr;
  /// Pointer to non-local old delayed incoming angular fluxes.
  double* nonlocal_delayed_incoming_psi_old = nullptr;
  /// Pointer to non-local outgoing offset.
  std::uint64_t* nonlocal_outgoing_offsets = nullptr;
  /// Pointer to non-local outgoing angular fluxes.
  double* nonlocal_outgoing_psi = nullptr;
  /// Pointer to boundary angular fluxes.
  double* boundary_psi = nullptr;
  /// Stride size (number of groups times number of angles).
  std::uint64_t stride_size;

  /// Get pointer to the incoming angular flux (if the face is not incoming, a nullptr is returned).
  constexpr double* GetIncomingFluxPointer(const AAHD_NodeIndex& node_index) const noexcept
  {
    // undefined case (parallel face) : all tests are true
    if (node_index.IsUndefined())
    {
      return nullptr;
    }
    // outgoing case : nullptr
    if (node_index.IsOutGoing())
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
        return nonlocal_delayed_incoming_psi_old +
               nonlocal_delayed_incoming_offsets[node_index.GetRank()] +
               node_index.GetIndex() * stride_size;
      }
      else
      {
        return nonlocal_incoming_psi + nonlocal_incoming_offsets[node_index.GetRank()] +
               node_index.GetIndex() * stride_size;
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
    if (!node_index.IsOutGoing())
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
      return nonlocal_outgoing_psi + nonlocal_outgoing_offsets[node_index.GetRank()] +
             node_index.GetIndex() * stride_size;
    }
  }

  /// Get pointer to the angular flux corresponding to the given node index.
  constexpr double* GetPointer(const AAHD_NodeIndex& node_index) const noexcept
  {
    // undefined case (parallel face) : all tests are true
    if (node_index.IsUndefined())
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
        if (node_index.IsOutGoing())
        {
          return delayed_local_psi + node_index.GetIndex() * stride_size;
        }
        else
        {
          return delayed_local_psi_old + node_index.GetIndex() * stride_size;
        }
      }
      else
      {
        return local_psi + node_index.GetIndex() * stride_size;
      }
    }
    // non-local case
    else
    {
      if (node_index.IsOutGoing())
      {
        return nonlocal_outgoing_psi + nonlocal_outgoing_offsets[node_index.GetRank()] +
               node_index.GetIndex() * stride_size;
      }
      else
      {
        if (node_index.IsDelayed())
        {
          return nonlocal_delayed_incoming_psi_old +
                 nonlocal_delayed_incoming_offsets[node_index.GetRank()] +
                 node_index.GetIndex() * stride_size;
        }
        else
        {
          return nonlocal_incoming_psi + nonlocal_incoming_offsets[node_index.GetRank()] +
                 node_index.GetIndex() * stride_size;
        }
      }
    }
  }
};

} // namespace opensn

/// Print face node.
std::ostream& operator<<(std::ostream& os, const opensn::AAHD_FaceNode& n);

/// Print non-local face node.
std::ostream& operator<<(std::ostream& os, const opensn::AAHD_NonLocalFaceNode& n);

/// Print node index.
std::ostream& operator<<(std::ostream& os, const opensn::AAHD_NodeIndex& e);
