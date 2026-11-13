// SPDX-FileCopyrightText: 2026 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep/communicators/cbc_message_parser.h"
#include <cstring>
#include <stdexcept>
#include <utility>

namespace opensn::detail
{
CBCMessageParser::CBCMessageParser(std::vector<Slot> normal,
                                   std::vector<Slot> delayed,
                                   std::span<const int> delayed_sources,
                                   int num_sources)
  : normal_(std::move(normal)), delayed_(std::move(delayed)), num_sources_(num_sources)
{
  if (num_sources < 0)
    throw std::logic_error("CBC parser has a negative source count");
  if (delayed_sources.size() > static_cast<std::size_t>(num_sources))
    throw std::logic_error("CBC parser has too many delayed sources");
  peers_.resize(delayed_sources.size());
  if (not delayed_sources.empty())
    delayed_peer_by_source_.assign(static_cast<std::size_t>(num_sources), peers_.size());
  for (std::size_t i = 0; i < delayed_sources.size(); ++i)
  {
    const int source = delayed_sources[i];
    if (source < 0 or source >= num_sources or delayed_peer_by_source_[source] != peers_.size())
      throw std::logic_error("CBC parser has an invalid/duplicate delayed source");
    delayed_peer_by_source_[source] = i;
  }
  for (const auto* slots : {&normal_, &delayed_})
    for (const auto& slot : *slots)
      if (slot.source < 0 or slot.source >= num_sources or slot.size == 0)
        throw std::logic_error("CBC parser has invalid face metadata");
  for (const auto& slot : delayed_)
    if (delayed_peer_by_source_.empty() or delayed_peer_by_source_[slot.source] == peers_.size())
      throw std::logic_error("CBC delayed face has no expected source");
  Reset();
}

void
CBCMessageParser::Reset()
{
  remaining_peers_ = peers_.size();
  for (auto& peer : peers_)
  {
    peer.complete = false;
    peer.remaining_faces = 0;
  }
  for (auto& slot : normal_)
    slot.received = 0;
  for (auto& slot : delayed_)
  {
    slot.received = 0;
    ++peers_[delayed_peer_by_source_[slot.source]].remaining_faces;
  }
}

void
CBCMessageParser::Parse(std::span<const char> packet,
                        int source,
                        bool delayed_stream,
                        const std::function<void(const Record&)>& store)
{
  if (source < 0 or source >= num_sources_ or packet.empty())
    throw std::logic_error("CBC packet has an invalid source or empty payload");
  const auto Read = [&]<typename T>()
  {
    if (packet.size() < sizeof(T))
      throw std::logic_error("Truncated CBC packet header");
    T value;
    std::memcpy(&value, packet.data(), sizeof(T));
    packet = packet.subspan(sizeof(T));
    return value;
  };
  while (not packet.empty())
  {
    const auto kind = Read.template operator()<std::uint8_t>();
    if (kind > static_cast<std::uint8_t>(CBCMessageKind::DELAYED_COMPLETION))
      throw std::logic_error("CBC packet has an invalid record kind");
    const auto peer_index = delayed_stream and not delayed_peer_by_source_.empty()
                              ? delayed_peer_by_source_[source]
                              : peers_.size();
    if (kind == static_cast<std::uint8_t>(CBCMessageKind::DELAYED_COMPLETION))
    {
      if (not delayed_stream or peer_index == peers_.size() or peers_[peer_index].complete or
          peers_[peer_index].remaining_faces != 0)
        throw std::logic_error("CBC delayed completion is unexpected, duplicate or premature");
      peers_[peer_index].complete = true;
      --remaining_peers_;
      continue;
    }
    const bool delayed = kind == static_cast<std::uint8_t>(CBCMessageKind::DELAYED_FACE_PSI) or
                         kind == static_cast<std::uint8_t>(CBCMessageKind::DELAYED_FACE_PSI_CHUNK);
    if (delayed != delayed_stream or
        (delayed and (peer_index == peers_.size() or peers_[peer_index].complete)))
      throw std::logic_error("CBC face record belongs to the wrong stream/source state");
    const auto index = Read.template operator()<std::size_t>();
    auto& slots = delayed ? delayed_ : normal_;
    if (index >= slots.size() or slots[index].source != source)
      throw std::logic_error("CBC face record has an invalid slot/source");
    auto& slot = slots[index];
    std::size_t total = slot.size, offset = 0, count = total;
    if (kind == static_cast<std::uint8_t>(CBCMessageKind::NORMAL_FACE_PSI_CHUNK) or
        kind == static_cast<std::uint8_t>(CBCMessageKind::DELAYED_FACE_PSI_CHUNK))
    {
      total = Read.template operator()<std::size_t>();
      offset = Read.template operator()<std::size_t>();
      count = Read.template operator()<std::size_t>();
    }
    // Division and subtraction avoid overflow even for hostile native-width length fields.
    if (total != slot.size or offset != slot.received or offset >= total or count == 0 or
        count > total - offset or count > packet.size() / sizeof(double))
      throw std::logic_error("CBC face record has invalid/truncated/duplicate fragment bounds");
    const bool complete = count == total - offset;
    store({delayed, index, total, offset, count, packet.data(), complete});
    slot.received += count;
    if (delayed and complete)
      --peers_[peer_index].remaining_faces;
    packet = packet.subspan(count * sizeof(double));
  }
}
} // namespace opensn::detail
