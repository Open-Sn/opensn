// SPDX-FileCopyrightText: 2026 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include <cstddef>
#include <cstdint>
#include <span>
#include <vector>

namespace opensn::detail
{
/// Collectively compare exact ordering/configuration descriptors, including empty descriptors.
/// Throws the same logic_error on every rank before any count-dependent CBC exchange.
void CheckCBCOrderingAgreement(std::span<const std::uint64_t> descriptor);

/// Validate a complete flattened CBC FAS stream before indexing any SPDS. Every ordinal occurs
/// exactly once; counts are even, endpoints are distinct communicator ranks, and pairs are unique.
std::vector<std::vector<int>>
ParseCBCSweepFAS(std::span<const int> payload, std::size_t num_orderings, int comm_size);
} // namespace opensn::detail
