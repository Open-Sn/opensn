// SPDX-FileCopyrightText: 2026 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include "hdf5.h"
#include <cstdint>
#include <functional>
#include <optional>
#include <string>
#include <vector>

namespace opensn
{

class DiscreteOrdinatesProblem;

struct UncollidedFluxData
{
  unsigned int num_groups = 0;
  unsigned int max_moment_order = 0;
  std::uint64_t global_cell_count = 0;
  double source_rate = 0.0;
  double outflow_rate = 0.0;
  std::vector<double> local_flux_moments;
};

class DiscreteOrdinatesProblemIO
{
public:
  static UncollidedFluxData ReadUncollidedFlux(const DiscreteOrdinatesProblem& do_problem,
                                               const std::string& file_name);

  static bool ReadRestartData(DiscreteOrdinatesProblem& do_problem,
                              hid_t file_id,
                              bool allow_transient_initialization_from_steady);

  static bool WriteRestartData(const DiscreteOrdinatesProblem& do_problem, hid_t file_id);

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
