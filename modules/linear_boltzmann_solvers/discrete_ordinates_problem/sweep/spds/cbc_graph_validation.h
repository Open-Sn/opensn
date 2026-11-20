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
/**
 * Label strongly connected components of the interpartition dependency graph collectively.
 * Only location zero receives the sparse graph. All ranks in mpi_comm participate.
 *
 * \param incoming_ranks Predecessor ranks, consistent with their outgoing neighbor lists.
 * \param outgoing_ranks Successor ranks, consistent with their incoming neighbor lists.
 * \return Identical component labels on every rank, with -1 for trivial components.
 * \throw std::logic_error If rank endpoints, counts or paired neighbor lists are invalid.
 */
std::vector<int> GetCBCCyclicRankComponents(const std::vector<int>& incoming_ranks,
                                            const std::vector<int>& outgoing_ranks);

/// Remote face endpoint used only during collective CBC no-lag cycle validation.
struct CBCBoundaryEdge
{
  /// Local cell index on this rank.
  std::size_t local_cell;
  /// Rank owning the adjacent cell.
  int remote_rank;
  /// Global ID of the adjacent cell on the remote rank.
  std::uint64_t remote_cell;
};

/**
 * Collectively check a directed cell graph whose owned subgraphs are DAGs.
 *
 * Every rank in mpi_comm participates and receives the same verdict. Ranks outside
 * nontrivial rank components may supply empty cell data. Exchange failures abort
 * because peers may be blocked.
 *
 * \param global_ids Owned cell global IDs, indexed by local cell ID.
 * \param successors Local successor IDs in an acyclic local graph, one row per owned cell.
 * \param incoming Remote incoming faces and their receiving local cells.
 * \param outgoing Remote outgoing faces and their sending local cells.
 * \param rank_components Identical rank-component labels on all ranks, or empty to compute them.
 * \return Identical empty string for an acyclic graph, or a contracted boundary-cycle witness.
 */
std::string CheckCBCGlobalCellGraph(const std::vector<std::uint64_t>& global_ids,
                                    const std::vector<std::vector<std::uint32_t>>& successors,
                                    const std::vector<CBCBoundaryEdge>& incoming,
                                    const std::vector<CBCBoundaryEdge>& outgoing,
                                    std::span<const int> rank_components = {});
} // namespace opensn::detail
