// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "physics/problems/linear_boltzmann/lbs_problem/groupset/lbs_groupset.h"
#include <vector>
#include <optional>

namespace opensn
{

class LBSProblem;

double ComputePointwisePhiChange(
  LBSProblem& lbs_problem,
  int groupset_id,
  std::optional<const std::reference_wrapper<std::vector<double>>> opt_phi_old = std::nullopt);

double ComputePointwisePhiChange(
  LBSProblem& lbs_problem,
  std::optional<const std::reference_wrapper<std::vector<double>>> opt_phi_old = std::nullopt);

double ComputePointwisePhiChange(
  LBSProblem& lbs_problem,
  std::vector<int> groupset_ids,
  std::optional<const std::reference_wrapper<std::vector<double>>> opt_phi_old = std::nullopt);

double ComputeL2PhiChange(
  LBSProblem& lbs_problem,
  std::optional<const std::reference_wrapper<std::vector<double>>> opt_phi_old = std::nullopt);

} // namespace opensn
