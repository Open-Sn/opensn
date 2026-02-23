// SPDX-FileCopyrightText: 2026 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include <functional>
#include <optional>
#include <string>
#include <vector>

namespace opensn
{

class DiscreteOrdinatesProblem;

class DiscreteOrdinatesProblemIO
{
public:
  /**
   * Write an angular flux vector to a file.
   *
   * \param do_problem Discrete ordinates problem
   * \param file_base File name base
   * \param opt_src Optional angular flux source vector
   */
  static void WriteAngularFluxes(
    DiscreteOrdinatesProblem& do_problem,
    const std::string& file_base,
    std::optional<const std::reference_wrapper<std::vector<std::vector<double>>>> opt_src =
      std::nullopt);

  /**
   * Read an angular flux vector from a file.
   *
   * \param do_problem Discrete ordinates problem
   * \param file_base File name base
   * \param opt_dest Optional angular flux destination vector
   */
  static void ReadAngularFluxes(
    DiscreteOrdinatesProblem& do_problem,
    const std::string& file_base,
    std::optional<std::reference_wrapper<std::vector<std::vector<double>>>> opt_dest =
      std::nullopt);
};

} // namespace opensn
