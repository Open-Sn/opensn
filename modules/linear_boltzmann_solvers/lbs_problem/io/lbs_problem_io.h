// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include <string>
#include <optional>
#include <vector>
#include <functional>

namespace opensn
{

class LBSProblem;
class DiscreteOrdinatesProblem;

class LBSSolverIO
{
public:
  /**
   * Write an angular flux vector to a file.
   *
   * \param lbs_problem LBS problem
   * \param file_base File name base
   * \param per_material Optional angular flux source vector
   */
  static void WriteAngularFluxes(
    DiscreteOrdinatesProblem& do_problem,
    const std::string& file_base,
    std::optional<const std::reference_wrapper<std::vector<std::vector<double>>>> opt_src =
      std::nullopt);

  /**
   * Read an angular flux vector from a file.
   *
   * \param lbs_problem LBS problem
   * \param file_base File name base
   * \param per_material Optional angular flux destination vector
   */
  static void ReadAngularFluxes(
    DiscreteOrdinatesProblem& do_problem,
    const std::string& file_base,
    std::optional<std::reference_wrapper<std::vector<std::vector<double>>>> opt_dest =
      std::nullopt);

  /**
   * Write a flux moments vector to a file.
   *
   * \param lbs_problem LBS problem
   * \param file_base File name base
   * \param per_material Optional flux moments source vector
   */
  static void WriteFluxMoments(
    LBSProblem& lbs_problem,
    const std::string& file_base,
    std::optional<const std::reference_wrapper<std::vector<double>>> opt_src = std::nullopt);

  /**
   * Read a flux moments vector from a file.
   *
   * \param lbs_problem LBS problem
   * \param file_base File name base
   * \param single_file Single data file or data file per rank?
   * \param per_material Optional flux moments destination vector
   */
  static void ReadFluxMoments(
    LBSProblem& lbs_problem,
    const std::string& file_base,
    bool single_file,
    std::optional<std::reference_wrapper<std::vector<double>>> opt_dest = std::nullopt);
};

} // namespace opensn
