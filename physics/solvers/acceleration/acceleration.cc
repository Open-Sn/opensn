// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "physics/solvers/acceleration/acceleration.h"
#include "physics/problems/linear_boltzmann/discrete_ordinates_problem/sweep/boundary/sweep_boundary.h"
#include "framework/materials/multi_group_xs/multi_group_xs.h"
#include "framework/logging/log_exceptions.h"
#include "framework/runtime.h"

namespace opensn
{

std::map<uint64_t, BoundaryCondition>
TranslateBCs(const std::map<uint64_t, std::shared_ptr<SweepBoundary>>& sweep_boundaries,
             bool vacuum_bcs_are_dirichlet)
{
  std::map<uint64_t, BoundaryCondition> bcs;
  for (auto& [bid, lbs_bndry] : sweep_boundaries)
  {
    if (lbs_bndry->GetType() == LBSBoundaryType::REFLECTING)
      bcs[bid] = {BCType::ROBIN, {0.0, 1.0, 0.0}};
    else if (lbs_bndry->GetType() == LBSBoundaryType::VACUUM)
      if (vacuum_bcs_are_dirichlet)
        bcs[bid] = {BCType::DIRICHLET, {0.0, 0.0, 0.0}};
      else
        bcs[bid] = {BCType::ROBIN, {0.25, 0.5}};
    else // dirichlet
      bcs[bid] = {BCType::DIRICHLET, {0.0, 0.0, 0.0}};
  }

  return bcs;
}

std::map<int, Multigroup_D_and_sigR>
PackGroupsetXS(const std::map<int, std::shared_ptr<MultiGroupXS>>& matid_to_xs_map,
               int first_grp_index,
               int last_group_index)
{
  const int num_gs_groups = last_group_index - first_grp_index + 1;
  OpenSnInvalidArgumentIf(num_gs_groups < 0, "last_grp_index must be >= first_grp_index");

  std::map<int, Multigroup_D_and_sigR> matid_2_mgxs_map;
  for (const auto& matid_xs_pair : matid_to_xs_map)
  {
    const auto& mat_id = matid_xs_pair.first;
    const auto& xs = matid_xs_pair.second;

    std::vector<double> D(num_gs_groups, 0.0);
    std::vector<double> sigma_r(num_gs_groups, 0.0);

    size_t g = 0;
    const auto& diffusion_coeff = xs->GetDiffusionCoefficient();
    const auto& sigma_removal = xs->GetSigmaRemoval();
    for (size_t gprime = first_grp_index; gprime <= last_group_index; ++gprime)
    {
      D[g] = diffusion_coeff[gprime];
      sigma_r[g] = sigma_removal[gprime];
      ++g;
    } // for g

    matid_2_mgxs_map.insert(std::make_pair(mat_id, Multigroup_D_and_sigR{D, sigma_r}));
  }

  return matid_2_mgxs_map;
}

} // namespace opensn
