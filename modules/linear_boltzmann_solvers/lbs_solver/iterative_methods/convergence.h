// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "modules/linear_boltzmann_solvers/lbs_solver/groupset/lbs_groupset.h"
#include <vector>
#include <optional>

namespace opensn
{

class LBSSolver;

double ComputePointwisePhiChange(
  LBSSolver& lbs_solver,
  int groupset_id,
  std::optional<const std::reference_wrapper<std::vector<double>>> opt_phi_old = std::nullopt);

double ComputePointwisePhiChange(
  LBSSolver& lbs_solver,
  std::optional<const std::reference_wrapper<std::vector<double>>> opt_phi_old = std::nullopt);

double ComputePointwisePhiChange(
  LBSSolver& lbs_solver,
  std::vector<int> groupset_ids,
  std::optional<const std::reference_wrapper<std::vector<double>>> opt_phi_old = std::nullopt);

double ComputeL2PhiChange(
  LBSSolver& lbs_solver,
  std::optional<const std::reference_wrapper<std::vector<double>>> opt_phi_old = std::nullopt);

} // namespace opensn