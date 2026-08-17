// SPDX-FileCopyrightText: 2026 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include <cstddef>
#include <cstdint>

namespace opensn
{

/// Message kinds shared by the host CBC serializer and communicator.
enum class CBC_MessageKind : std::uint8_t
{
  NORMAL_FACE_PSI = 0,        ///< Complete ordinary face values.
  DELAYED_FACE_PSI = 1,       ///< Complete delayed face values.
  NORMAL_FACE_PSI_CHUNK = 2,  ///< Fragment of ordinary face values.
  DELAYED_FACE_PSI_CHUNK = 3, ///< Fragment of delayed face values.
  DELAYED_COMPLETION = 4      ///< End of a location's delayed messages for this sweep.
};

/**
 * Receiving message details for one host CBC nonlocal face.
 * Face metadata is fixed for an angle set. Reset clears receive progress.
 *
 * All ranks must use the same native field representation. Fields are packed without
 * padding. A face record holds a uint8_t kind, a size_t receiver slot and double values.
 * Fragments insert three size_t fields before the values (total value count, offset and
 * fragment count). A completion record holds only its kind. Record kinds must match the
 * ordinary or delayed MPI tag.
 */
struct CBC_MessageDetails
{
  /// Expected source rank in the sweep communicator.
  int peer = -1;
  /// Expected number of double values in the complete face.
  std::size_t size = 0;
  /// Number of contiguous values stored in this sweep.
  std::size_t received = 0;
};

} // namespace opensn
