// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include "hdf5.h"
#include <cstddef>
#include <cstdint>
#include <string>
#include <optional>
#include <vector>
#include <functional>
#include <filesystem>

namespace opensn
{

class LBSProblem;

enum class FluxSolutionRole : unsigned int
{
  UNSPECIFIED = 0,
  FORWARD = 1,
  ADJOINT = 2
};

enum class AdjointFluxOrientation : unsigned int
{
  UNSPECIFIED = 0,
  SOLVER = 1,
  PHYSICAL = 2
};

class LBSSolverIO
{
public:
  static bool ReadRestartData(LBSProblem& lbs_problem,
                              const std::function<bool(hid_t)>& extra_reader = {},
                              const std::filesystem::path& read_path = {},
                              bool allow_transient_initialization_from_steady = false);

  static bool WriteRestartData(LBSProblem& lbs_problem,
                               const std::function<bool(hid_t)>& extra_writer = {});

  /**
   * Write a flux moments vector to one file per rank, `<file_base><rank>.h5`.
   *
   * Each file holds the rank's local cells: their global ids, node counts, and node locations in
   * the "mesh" group, and the flux moments in the dataset "values", ordered by cell, node, moment,
   * and group. The root attributes record the number of moments and groups, the solution role and
   * adjoint orientation, `num_partitions` (uint64), the number of files in the set, and
   * `file_set_id` (uint64), a random identifier shared by the files of this write.
   *
   * Collective: every rank must pass the same `file_base`, `role`, and `orientation`. Every rank
   * checks its source before any rank replaces its file, so a rejected source leaves an existing
   * set untouched; then every rank writes its file. A failure on any rank, such as a source that
   * does not match the problem or a failed write, raises std::runtime_error on every rank, since
   * the set is complete only if every file is. Inconsistent `role` and `orientation` raise
   * std::invalid_argument before any file is written.
   *
   * \param lbs_problem LBS problem
   * \param file_base File name base
   * \param opt_src Optional source vector to write instead of the problem-owned flux moments.
   * \param role Solution role. For problem-owned data, unspecified is inferred from the problem.
   * \param orientation Adjoint orientation. For problem-owned data, unspecified is inferred.
   */
  static void WriteFluxMoments(
    LBSProblem& lbs_problem,
    const std::string& file_base,
    std::optional<const std::reference_wrapper<std::vector<double>>> opt_src = std::nullopt,
    FluxSolutionRole role = FluxSolutionRole::UNSPECIFIED,
    AdjointFluxOrientation orientation = AdjointFluxOrientation::UNSPECIFIED);

  /**
   * Read a flux moments vector from a set of files written by WriteFluxMoments.
   *
   * \param lbs_problem LBS problem
   * \param file_base File name base
   * \param single_file Read the single file `<file_base>.h5`, which must not record a set of more
   * than one file, instead of the set `<file_base><k>.h5`.
   * \param opt_dest Optional destination vector to populate instead of the problem-owned flux
   * moments.
   * \param expected_role Required solution role, or unspecified to accept any role.
   * \param expected_orientation Required orientation of a file that records the adjoint role, or
   * unspecified to accept any. Files of other roles, and files without provenance, pass.
   *
   * The files may have been written by any number of ranks with any partitioning of the same mesh;
   * ReadCellwiseFileSet sends each cell's values to the rank that owns it, so every rank of this
   * run receives exactly its own cells. Numbered files without `num_partitions`, written by earlier
   * versions, form a set of as many files as this run has ranks; a single file (`single_file`) is
   * always one file and can be read on any number of ranks. Collective: a file that cannot be read
   * or does not match the problem, a local cell missing from the set, files of the set that record
   * different write identifiers, or files of the set that record different solution provenance
   * (roles, adjoint orientations, or whether they record them) raise std::runtime_error on every
   * rank, and the destination is left unchanged. Returns the adjoint orientation recorded in the
   * files, the same on every rank, unspecified for a forward file or a file without provenance.
   */
  static AdjointFluxOrientation ReadFluxMoments(
    LBSProblem& lbs_problem,
    const std::string& file_base,
    bool single_file,
    std::optional<std::reference_wrapper<std::vector<double>>> opt_dest = std::nullopt,
    FluxSolutionRole expected_role = FluxSolutionRole::UNSPECIFIED,
    AdjointFluxOrientation expected_orientation = AdjointFluxOrientation::UNSPECIFIED);
};

} // namespace opensn
