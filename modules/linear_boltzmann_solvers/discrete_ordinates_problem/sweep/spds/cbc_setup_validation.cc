// SPDX-FileCopyrightText: 2026 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep/spds/cbc_setup_validation.h"
#include "framework/runtime.h"
#include <algorithm>
#include <array>
#include <limits>
#include <set>
#include <stdexcept>
#include <utility>

namespace opensn::detail
{
void
CheckCBCOrderingAgreement(std::span<const std::uint64_t> descriptor)
{
  const auto n = descriptor.size(), max = std::numeric_limits<std::size_t>::max();
  const std::array<std::size_t, 2> local{n, max - n};
  std::array<std::size_t, 2> bounds{};
  mpi_comm.all_reduce(
    local.data(), static_cast<int>(local.size()), bounds.data(), mpi::op::max<std::size_t>());
  if (bounds[0] != max - bounds[1])
    throw std::logic_error("CBC SPDS count/configuration differs across MPI ranks");
  // No hash collision can hide a configuration mismatch. Broadcast in chunks bounded by the MPI
  // count type; the comparison storage is O(ordering descriptors), independent of cell counts.
  std::vector<std::uint64_t> reference(descriptor.begin(), descriptor.end());
  std::size_t offset = 0;
  while (offset < n)
  {
    const auto count = static_cast<int>(
      std::min(n - offset, static_cast<std::size_t>(std::numeric_limits<int>::max())));
    mpi_comm.broadcast(reference.data() + offset, count, 0);
    offset += static_cast<std::size_t>(count);
  }
  int differs = not std::equal(descriptor.begin(), descriptor.end(), reference.begin());
  int any_differs = 0;
  mpi_comm.all_reduce(differs, any_differs, mpi::op::max<int>());
  if (any_differs)
    throw std::logic_error("CBC SPDS ordering/configuration differs across MPI ranks");
}

std::vector<std::vector<int>>
ParseCBCSweepFAS(std::span<const int> payload, std::size_t num_orderings, int comm_size)
{
  std::vector<std::vector<int>> result(num_orderings);
  std::vector<bool> received(num_orderings, false);
  std::size_t offset = 0;
  while (offset < payload.size())
  {
    if (payload.size() - offset < 2)
      throw std::logic_error("Truncated CBC sweep FAS header");
    const int ordinal = payload[offset++], count = payload[offset++];
    if (ordinal < 0 or std::cmp_greater_equal(ordinal, num_orderings) or count < 0 or
        count % 2 != 0 or static_cast<std::size_t>(count) > payload.size() - offset)
      throw std::logic_error("Invalid CBC sweep FAS ordinal/count");
    if (received[ordinal])
      throw std::logic_error("Duplicate CBC sweep FAS ordinal");
    received[ordinal] = true;
    std::set<std::pair<int, int>> pairs;
    for (int i = 0; i < count; i += 2)
    {
      const int source = payload[offset + i], target = payload[offset + i + 1];
      if (source < 0 or target < 0 or source >= comm_size or target >= comm_size or
          source == target or not pairs.emplace(source, target).second)
        throw std::logic_error("Invalid/duplicate CBC sweep FAS rank pair");
    }
    result[ordinal].assign(payload.data() + offset, payload.data() + offset + count);
    offset += static_cast<std::size_t>(count);
  }
  if (std::find(received.begin(), received.end(), false) != received.end())
    throw std::logic_error("Missing CBC sweep FAS ordinal");
  return result;
}
} // namespace opensn::detail
