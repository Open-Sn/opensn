// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include <string>
#include <optional>
#include <vector>
#include <functional>

namespace opensn
{

class LBSSolver;

class LBSSolverIO
{
public:
  /**
   * Write an angular flux vector to a file.
   *
   * \param lbs_solver LBS solver
   * \param file_base File name stem
   * \param per_material Optional angular flux source vector
   */
  static void WriteAngularFluxes(
    LBSSolver& lbs_solver,
    const std::string& file_stem,
    std::optional<const std::reference_wrapper<std::vector<std::vector<double>>>> opt_src =
      std::nullopt);

  /**
   * Read an angular flux vector from a file.
   *
   * \param lbs_solver LBS solver
   * \param file_base File name stem
   * \param per_material Optional angular flux destination vector
   */
  static void ReadAngularFluxes(
    LBSSolver& lbs_solver,
    const std::string& file_stem,
    std::optional<std::reference_wrapper<std::vector<std::vector<double>>>> opt_dest =
      std::nullopt);

  /**
   * Write an angular flux groupset vector to a file.
   *
   * \param lbs_solver LBS solver
   * \param groupset_id Energy groupset id
   * \param file_base File name stem
   * \param per_material Optional angular flux source vector
   */
  static void WriteGroupsetAngularFluxes(
    LBSSolver& lbs_solver,
    const int groupset_id,
    const std::string& file_stem,
    std::optional<const std::reference_wrapper<std::vector<double>>> opt_src = std::nullopt);

  /**
   * Read an angular flux groupset vector from a file.
   *
   * \param lbs_solver LBS solver
   * \param groupset_id Energy groupset id
   * \param file_base File name stem
   * \param per_material Optional angular flux destination vector
   */
  static void ReadGroupsetAngularFluxes(
    LBSSolver& lbs_solver,
    const int groupset_id,
    const std::string& file_base,
    std::optional<std::reference_wrapper<std::vector<double>>> opt_dest = std::nullopt);

  /**
   * Write a flux moments vector to a file.
   *
   * \param lbs_solver LBS solver
   * \param file_base File name stem
   * \param per_material Optional flux moments source vector
   */
  static void WriteFluxMoments(
    LBSSolver& lbs_solver,
    const std::string& file_stem,
    std::optional<const std::reference_wrapper<std::vector<double>>> opt_src = std::nullopt);

  /**
   * Read a flux moments vector from a file.
   *
   * \param lbs_solver LBS solver
   * \param file_base File name stem
   * \param single_file Single data file or data file per rank?
   * \param per_material Optional flux moments destination vector
   */
  static void ReadFluxMoments(
    LBSSolver& lbs_solver,
    const std::string& file_stem,
    bool single_file,
    std::optional<std::reference_wrapper<std::vector<double>>> opt_dest = std::nullopt);
};

} // namespace opensn
