// SPDX-FileCopyrightText: 2026 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include <cstddef>
#include <cstdint>
#include <functional>
#include <span>
#include <vector>

namespace opensn::detail
{
/// Native wire discriminants shared by the CBC serializer and checked parser.
enum class CBCMessageKind : std::uint8_t
{
  NORMAL_FACE_PSI = 0,
  DELAYED_FACE_PSI = 1,
  NORMAL_FACE_PSI_CHUNK = 2,
  DELAYED_FACE_PSI_CHUNK = 3,
  DELAYED_COMPLETION = 4
};

/// Validated receiver of the host CBC native-width, same-executable MPI wire protocol.
class CBCMessageParser
{
public:
  struct Slot
  {
    int source = -1;
    std::size_t size = 0;
    std::size_t received = 0;
  };
  struct Record
  {
    bool delayed;
    std::size_t slot, total, offset, count;
    /// Borrowed packet storage, valid only during the synchronous store callback.
    const char* payload;
    bool complete;
  };

  /// Metadata is fixed for one angle set. Sources are ranks in the sweep communicator.
  CBCMessageParser(std::vector<Slot> normal,
                   std::vector<Slot> delayed,
                   std::span<const int> delayed_sources,
                   int num_sources);
  /// Reset progress only after draining both streams; fixed slot metadata survives.
  void Reset();
  /// Parse a packet and call store only for validated records. Throws logic_error on malformed
  /// input. Records from each source/tag must preserve sender order (the MPI ordering contract).
  /// Fragment offsets must equal the received prefix, detecting overlaps/duplicates in O(1).
  void Parse(std::span<const char> packet,
             int source,
             bool delayed_stream,
             const std::function<void(const Record&)>& store);
  bool DelayedComplete() const noexcept { return remaining_peers_ == 0; }

private:
  struct Peer
  {
    bool complete = false;
    std::size_t remaining_faces = 0;
  };
  std::vector<Slot> normal_, delayed_;
  std::vector<Peer> peers_;
  /// As in the original communicator, allocate a rank-index table only with delayed sources.
  std::vector<std::size_t> delayed_peer_by_source_;
  int num_sources_ = 0;
  std::size_t remaining_peers_ = 0;
};
} // namespace opensn::detail
