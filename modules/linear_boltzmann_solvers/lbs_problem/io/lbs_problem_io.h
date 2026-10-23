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
class Mesh;
class SpatialDiscretization;

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
   * Write a flux moments vector to a file.
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
   * Read a flux moments vector from a file.
   *
   * \param lbs_problem LBS problem
   * \param file_base File name base
   * \param single_file Single data file or data file per rank?
   * \param opt_dest Optional destination vector to populate instead of the problem-owned flux
   * moments.
   * \param expected_role Required solution role, or unspecified to accept any role.
   * \param expected_orientation Required adjoint orientation, or unspecified to accept any.
   *
   * The file must have been written with the same rank count and partition. Collective: a file
   * that cannot be read or does not match the problem raises an error on every rank, and the
   * destination is left unchanged.
   */
  static void ReadFluxMoments(
    LBSProblem& lbs_problem,
    const std::string& file_base,
    bool single_file,
    std::optional<std::reference_wrapper<std::vector<double>>> opt_dest = std::nullopt,
    FluxSolutionRole expected_role = FluxSolutionRole::UNSPECIFIED,
    AdjointFluxOrientation expected_orientation = AdjointFluxOrientation::UNSPECIFIED);

  /// The cells of a flux-moment or angular-flux file and the local node of each file node.
  struct FileMeshNodes
  {
    std::vector<std::uint64_t> cell_ids;
    /// For each file cell, the local cell node matching each file node, in file order.
    std::vector<std::vector<std::size_t>> node_maps;
  };

  /**
   * Reads and validates the "mesh" group of a flux-moment or angular-flux file. Every file cell
   * must be a distinct local cell whose nodes match the file node locations one to one, and the
   * dataset sizes must match the metadata. Throws on the calling rank only.
   */
  static FileMeshNodes ReadFileMeshNodes(hid_t file,
                                         const std::string& file_name,
                                         const Mesh& grid,
                                         const SpatialDiscretization& discretization);
};

} // namespace opensn
