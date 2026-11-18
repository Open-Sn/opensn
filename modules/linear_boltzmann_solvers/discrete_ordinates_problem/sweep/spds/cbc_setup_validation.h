// SPDX-FileCopyrightText: 2026 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include <cstddef>
#include <cstdint>
#include <span>
#include <vector>

namespace opensn::detail
{
/**
 * Compare CBC setup descriptions collectively on mpi_comm, including empty descriptions.
 *
 * \param descriptor Serialized ordering IDs, cycle settings and direction data for this rank.
 * \throw std::logic_error On every rank if any descriptions differ.
 */
void CheckCBCOrderingAgreement(std::span<const std::uint64_t> descriptor);

/**
 * Validate a complete serialized CBC feedback arc set before indexing any SPDS.
 * Every ordering occurs exactly once; counts are even and rank pairs are distinct.
 *
 * \param payload Flattened ordering index, endpoint count and endpoint pairs for each ordering.
 * \param num_orderings Number of sweep orderings expected in the payload.
 * \param comm_size Number of ranks defining the valid endpoint range.
 * \return Removed rank pairs grouped by ordering index.
 * \throw std::logic_error If an ordering, count or endpoint pair is missing or invalid.
 */
std::vector<std::vector<int>>
ParseCBCSweepFAS(std::span<const int> payload, std::size_t num_orderings, int comm_size);
} // namespace opensn::detail
