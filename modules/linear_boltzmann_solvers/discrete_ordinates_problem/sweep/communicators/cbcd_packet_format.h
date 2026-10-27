// SPDX-FileCopyrightText: 2026 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include <cstddef>
#include <cstdint>
#include <limits>
#include <stdexcept>

namespace opensn::cbcd_packet
{

/**
 * \name Device CBC packet format
 * A packet sends the angular fluxes of one or more outgoing nonlocal faces to one destination.
 * It is sent as two messages:
 * - a host control message of words (a header followed by one record per face)
 * - a device payload of the faces' angular-flux values (in record order)
 * The control message is sent first and on its own tag so that the receiver can post the payload
 * receive and prepare the faces' addresses while the payload is still being packed and sent.
 * \{
 */

/// Header word: sequence number of the packet among the packets from the sender to the receiver.
constexpr std::size_t HEADER_SEQUENCE = 0;
/// Header word: number of payload values.
constexpr std::size_t HEADER_NUM_VALUES = 1;
/// Header word: number of face records.
constexpr std::size_t HEADER_NUM_FACES = 2;
/// Number of header words.
constexpr std::size_t HEADER_WORDS = 3;

/// Record word: index of the sending angle set.
constexpr std::size_t RECORD_ANGLE_SET = 0;
/// Record word: index of the face among the faces the receiver gets from that angle set.
constexpr std::size_t RECORD_DESTINATION_FACE = 1;
/// Record word: number of values of the face.
constexpr std::size_t RECORD_NUM_VALUES = 2;
/// Record word: offset of the face's first value in the payload.
constexpr std::size_t RECORD_PAYLOAD_OFFSET = 3;
/// Number of words per face record.
constexpr std::size_t RECORD_WORDS = 4;

/// Return the number of control words of a packet with `num_faces` faces.
constexpr std::size_t
ControlWords(std::size_t num_faces)
{
  return HEADER_WORDS + RECORD_WORDS * num_faces;
}

/// Largest element count of one MPI message.
constexpr std::size_t MPI_MAX_COUNT = static_cast<std::size_t>(std::numeric_limits<int>::max());

/// Largest number of faces per packet so that the control message fits MPI_MAX_COUNT words.
constexpr std::size_t MAX_FACES = (MPI_MAX_COUNT - HEADER_WORDS) / RECORD_WORDS;

/**
 * Upper bound of the control words of all packets carrying `num_faces` faces in total.
 * Every packet carries at least one face, so there are at most `num_faces` headers.
 */
constexpr std::size_t
MaxControlWords(std::size_t num_faces)
{
  return (HEADER_WORDS + RECORD_WORDS) * num_faces;
}

/// Message-tag channels of a groupset.
constexpr std::size_t TAG_STRIDE = 2;
/// Tag channel of the control messages.
constexpr std::size_t CONTROL_TAG_OFFSET = 0;
/// Tag channel of the payload messages.
constexpr std::size_t PAYLOAD_TAG_OFFSET = 1;
/** \} */

/// Return `left + right`.
inline std::size_t
CheckedAdd(std::size_t left, std::size_t right)
{
  if (right > std::numeric_limits<std::size_t>::max() - left)
    throw std::overflow_error("CBCD device-packet storage exceeds size_t");
  return left + right;
}

/// Return `left * right`.
inline std::size_t
CheckedMultiply(std::size_t left, std::size_t right)
{
  if (right != 0 and left > std::numeric_limits<std::size_t>::max() / right)
    throw std::overflow_error("CBCD device-packet storage exceeds size_t");
  return left * right;
}

} // namespace opensn::cbcd_packet
