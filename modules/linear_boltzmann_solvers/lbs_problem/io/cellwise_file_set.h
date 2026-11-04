// SPDX-FileCopyrightText: 2026 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include "modules/linear_boltzmann_solvers/lbs_problem/io/lbs_problem_io.h"
#include "framework/runtime.h"
#include "hdf5.h"
#include <cstddef>
#include <cstdint>
#include <functional>
#include <optional>
#include <string>
#include <vector>

namespace opensn
{

class Cell;
class Mesh;
class SpatialDiscretization;

/// A dataset of a cell-wise file holding `stride` values per file node, ordered by cell and node.
struct CellwiseDataset
{
  std::string name;
  std::uint64_t stride = 0;
};

/// Default bound on the packed bytes each rank sends per exchange while reading a cell-wise set.
const std::uint64_t cellwise_redistribution_bytes = std::uint64_t{64} << 20;

/**
 * Return a random identifier for the files of one write, the same on every rank of `comm`, drawn
 * on rank 0. Collective. Each file of a set records it, so that a reader rejects a set that
 * combines files of different writes. Throws std::runtime_error on every rank if rank 0 cannot
 * draw it.
 */
std::uint64_t MakeFileSetId(const mpi::Communicator& comm);

/**
 * Write the attributes `num_partitions`, the size of `comm`, and `file_set_id` (uint64), the
 * identifier of this write from MakeFileSetId, and the "mesh" group that ReadCellwiseFileSet reads:
 * the global ids, node counts, and node locations of the local cells of `grid`, in local cell
 * order. A cell-wise dataset written after it stores the values of each local node in the same
 * order. Throws std::logic_error on the calling rank if a write fails.
 */
void WriteCellwiseFileMesh(hid_t file,
                           const Mesh& grid,
                           const SpatialDiscretization& discretization,
                           const mpi::Communicator& comm,
                           std::uint64_t file_set_id);

/**
 * The solution provenance that a flux-moment or angular-flux file records: whether it has the
 * provenance attributes, and their solution role and adjoint orientation. Both are unspecified for
 * a file without the attributes, and the orientation is unspecified for a non-adjoint file. A file
 * may also record an unspecified role, as one written from a supplied vector does, so `tagged`
 * tells it apart from a file without the attributes.
 */
struct SolutionProvenance
{
  bool tagged = false;
  FluxSolutionRole role = FluxSolutionRole::UNSPECIFIED;
  AdjointFluxOrientation orientation = AdjointFluxOrientation::UNSPECIFIED;

  bool operator==(const SolutionProvenance& other) const = default;
};

/**
 * Check the solution provenance of a flux-moment or angular-flux file: the `solution_role` and
 * `adjoint_orientation` attributes must be both present or both absent, valid, and consistent; they
 * must match `expected_role` unless it is unspecified; an adjoint file must match
 * `expected_orientation` unless it is unspecified, while files of other roles and files without
 * provenance pass; and when the file replaces the state of a problem, `problem_is_adjoint` is the
 * problem's mode, a recorded role must match it, and an adjoint must be in the solver orientation.
 * Returns the recorded provenance. Throws std::logic_error on the calling rank.
 */
SolutionProvenance CheckSolutionProvenance(hid_t file,
                                           const std::string& file_name,
                                           FluxSolutionRole expected_role,
                                           AdjointFluxOrientation expected_orientation,
                                           std::optional<bool> problem_is_adjoint);

/**
 * Checks that every file of a cell-wise set records the same solution provenance, so that a set
 * never combines, for example, forward and adjoint files, solver- and physically oriented adjoint
 * files, or files with and without the provenance attributes, including files that record an
 * unspecified role.
 */
class SetProvenanceCheck
{
public:
  /// Add the provenance of a file this rank read. Throws std::logic_error on the calling rank if
  /// it differs from that of an earlier file of this rank.
  void Add(const SolutionProvenance& provenance, const std::string& file_name);

  /**
   * Return the provenance of the set. Collective over `comm`; a rank that read no file takes it
   * from the others. Throws std::runtime_error on every rank, with `operation` in its message, if
   * the files of the set record different provenance.
   */
  SolutionProvenance Agree(const mpi::Communicator& comm, const std::string& operation) const;

private:
  std::optional<SolutionProvenance> local_;
};

/**
 * Read a set of cell-wise files, such as those written by LBSSolverIO::WriteFluxMoments and
 * DiscreteOrdinatesProblemIO::WriteAngularFluxes, onto the cells of `grid` that this rank owns.
 *
 * The set is `<file_base>0.h5`, ..., `<file_base>(N-1).h5`, or the single file `<file_base>.h5` if
 * `single_file` is true. Each file stores the global ids, node counts, and node locations of its
 * cells in its "mesh" group, and each dataset of `datasets` stores `stride` values per node in
 * that order. N is the `num_partitions` attribute of the files, at most INT_MAX. Numbered files
 * written before it was recorded form a set of as many files as `comm` has ranks, so they must be
 * read on the rank count that wrote them. A single file is always one file: with or without
 * `num_partitions`, which must then be 1, it can be read on any number of ranks.
 *
 * The files may have been written by any number of ranks with any partitioning of the same mesh,
 * but must come from one write: every file that records `num_partitions` also records the
 * `file_set_id` of its write (see MakeFileSetId), and a set whose files record different ids is
 * rejected. Older files without `num_partitions` carry no id, so a set of them is not checked. A
 * file may hold at most the cells and nodes of the whole mesh; larger counts are rejected before
 * anything is allocated.
 *
 * File k is read by rank k mod P, one file at a time, and each cell's node locations and values
 * are sent to the rank that owns the cell, found with an OwnerDirectory, in exchanges in which each
 * rank sends at most `max_send_bytes` packed bytes. A cell larger than that is sent alone, so the
 * bound is not a hard memory limit. A rank may receive up to P times the bound, so the bound used
 * is min(`max_send_bytes`, 8 * floor(INT_MAX / P)) bytes (with the default, the second term is the
 * smaller on more than 256 ranks), which keeps every rank's receipt in one exchange within the MPI
 * int count range when every cell fits the bound. Cells
 * larger than the bound, sent together to one rank, can still exceed that range; the read then
 * fails on every rank. Before its streams, the owners of a file's cells are looked up in one
 * exchange outside the bound, of 24 bytes per file cell each way. A rank holds one file, one of its
 * datasets, and one exchange at a time, in addition to its result.
 *
 * Every rank calls `validate_file` on each file it opens, before its data are read, and calls
 * `scatter(dataset, cell, node, values)` for every node of every local cell: `values` are the
 * `stride` values of dataset `dataset` at local node `node` of `cell`. Both callbacks are
 * rank-local and run a different number of times on each rank, so they must not make MPI calls;
 * an exception they throw is raised on every rank. File nodes are matched to local nodes by
 * location, so their order within a cell may differ: first at exactly the same location, then
 * within a small tolerance when that match is unique.
 *
 * `max_send_bytes` must be positive; zero throws std::invalid_argument before any file is opened.
 *
 * Collective over `comm`. Throws std::runtime_error on every rank, with `operation` and the lowest
 * failing rank in its message, if a file is missing, malformed, inconsistent with the set, contains
 * non-finite values, fails `validate_file`, or comes from another write than file 0; if a file
 * holds more cells or nodes than the mesh; if a file cell is not a cell of `grid`, appears in
 * more than one file, or has node locations that differ from the mesh or cannot be matched
 * unambiguously; if a local cell is in no file; if one exchange would exceed the MPI int count
 * range; or if a rank-local step, such as an allocation, fails. A failure may leave `scatter`
 * called for some nodes, so callers scatter into a temporary.
 */
void ReadCellwiseFileSet(
  const std::string& file_base,
  bool single_file,
  const std::vector<CellwiseDataset>& datasets,
  const Mesh& grid,
  const SpatialDiscretization& discretization,
  const std::function<void(hid_t file, const std::string& file_name)>& validate_file,
  const std::function<
    void(std::size_t dataset, const Cell& cell, std::size_t node, const double* values)>& scatter,
  const mpi::Communicator& comm,
  const std::string& operation,
  std::uint64_t max_send_bytes = cellwise_redistribution_bytes);

} // namespace opensn
