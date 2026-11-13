// SPDX-FileCopyrightText: 2026 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include <cstddef>
#include <cstdint>
#include <string>
#include <span>
#include <vector>

namespace opensn::detail
{
/// Collectively return identical rank SCC labels, with -1 for trivial components. Only rank zero
/// receives the sparse quotient graph; incoming/outgoing neighbor lists must agree globally.
std::vector<int> GetCBCCyclicRankComponents(const std::vector<int>& incoming_ranks,
                                            const std::vector<int>& outgoing_ranks);

/// Remote face endpoint used only during collective CBC no-lag cycle validation.
struct CBCBoundaryEdge
{
  std::size_t local_cell;
  int remote_rank;
  std::uint64_t remote_cell;
};

/**
 * Collectively check a directed cell graph whose owned subgraphs are DAGs.
 *
 * Every rank participates, including empty ranks. Local IDs index global_ids and successors.
 * Incoming/outgoing contain the remote faces with their receiving/sending local cell respectively.
 * Returns the same empty string (acyclic) or boundary-cycle witness on all ranks. Uses mpi_comm.
 * rank_components may hold the labels returned by GetCBCCyclicRankComponents, identically on
 * every rank. With an empty span they are computed from the boundary edges. Ranks outside
 * nontrivial components may supply empty cell data. Exchange failures abort before communicator
 * teardown, since peers may still be blocked. No cell graph is gathered. See groupsets.rst for the
 * elimination proof, witness contraction and setup-only costs. Each cell/local edge is visited
 * once; remote faces are validated and released without forming transitive reachability pairs.
 */
std::string CheckCBCGlobalCellGraph(const std::vector<std::uint64_t>& global_ids,
                                    const std::vector<std::vector<std::uint32_t>>& successors,
                                    const std::vector<CBCBoundaryEdge>& incoming,
                                    const std::vector<CBCBoundaryEdge>& outgoing,
                                    std::span<const int> rank_components = {});
} // namespace opensn::detail
