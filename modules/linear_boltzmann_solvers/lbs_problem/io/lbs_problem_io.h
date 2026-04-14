// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include "hdf5.h"
#include <string>
#include <optional>
#include <vector>
#include <functional>

namespace opensn
{

class LBSProblem;

class LBSSolverIO
{
public:
  static bool ReadRestartData(LBSProblem& lbs_problem,
                              const std::function<bool(hid_t)>& extra_reader = {});

  static bool WriteRestartData(LBSProblem& lbs_problem,
                               const std::function<bool(hid_t)>& extra_writer = {});

  /**
   * Write a flux moments vector to a file.
   *
   * \param lbs_problem LBS problem
   * \param file_base File name base
   * \param opt_src Optional source vector to write instead of the problem-owned flux moments.
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
   * \param opt_dest Optional destination vector to populate instead of the problem-owned flux
   * moments.
   */
  static void ReadFluxMoments(
    LBSProblem& lbs_problem,
    const std::string& file_base,
    bool single_file,
    std::optional<std::reference_wrapper<std::vector<double>>> opt_dest = std::nullopt);
};

} // namespace opensn
