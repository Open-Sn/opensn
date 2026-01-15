// SPDX-FileCopyrightText: 2025 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep/fluds/device_structs.h"

namespace opensn
{

using CBCD_FaceNode = DeviceFaceNode;

using CBCD_NodeIndex = DeviceNodeIndex;

/**
 * \brief Set of device pointers to local, boundary, and non-local buffers for CBCD FLUDS.
 */
struct CBCD_FLUDSPointerSet
{
  /// Pointer to device local angular flux buffer.
  double* local_psi = nullptr;

  /// Pointer to device incoming boundary angular flux buffer.
  double* incoming_boundary_psi = nullptr;

  /// Pointer to device outgoing boundary angular flux buffer.
  double* outgoing_boundary_psi = nullptr;

  /// Pointer to device incoming non-local angular flux buffer.
  double* incoming_nonlocal_psi = nullptr;

  /// Pointer to device outgoing non-local angular flux buffer.
  double* outgoing_nonlocal_psi = nullptr;

  // Stride size (number of groups times number of angles).
  std::uint64_t stride_size = 0;

  /// Get pointer to the incoming angular flux (if the face is not incoming, a nullptr is returned).
  constexpr double* GetIncomingFluxPointer(const CBCD_NodeIndex& node_index) const noexcept
  {
    // Undefined case (corresponds to a parallel face)
    if (node_index.IsUndefined())
    {
      return nullptr;
    }
    // Outgoing case : nullptr
    if (node_index.IsOutgoing())
    {
      return nullptr;
    }
    // Incoming boundary case
    if (node_index.IsBoundary())
    {
      return incoming_boundary_psi + node_index.GetIndex() * stride_size;
    }
    // Incoming local case
    if (node_index.IsLocal())
    {
      return local_psi + node_index.GetIndex() * stride_size;
    }
    // Incoming non-local case
    else
    {
      return incoming_nonlocal_psi + node_index.GetIndex() * stride_size;
    }
  }

  /// Get pointer to the outgoing angular flux (if the face is not outgoing, a nullptr is returned).
  constexpr double* GetOutgoingFluxPointer(const CBCD_NodeIndex& node_index) const noexcept
  {
    // Undefined case (corresponds to a parallel face)
    if (node_index.IsUndefined())
    {
      return nullptr;
    }
    // Incoming case : nullptr
    if (!node_index.IsOutgoing())
    {
      return nullptr;
    }
    // Outgoing boundary case
    if (node_index.IsBoundary())
    {
      return outgoing_boundary_psi + node_index.GetIndex() * stride_size;
    }
    // Outgoing local case
    if (node_index.IsLocal())
    {
      return local_psi + node_index.GetIndex() * stride_size;
    }
    // Outgoing non-local case
    else
    {
      return outgoing_nonlocal_psi + node_index.GetIndex() * stride_size;
    }
  }
};

} // namespace opensn