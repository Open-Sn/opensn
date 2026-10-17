// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/acceleration/acceleration.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep/boundary/sweep_boundary.h"
#include "framework/materials/multi_group_xs/multi_group_xs.h"
#include "framework/utils/error.h"
#include "framework/runtime.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/discrete_ordinates_problem.h"

namespace opensn
{

std::map<uint64_t, BoundaryCondition>
TranslateBCs(const std::map<uint64_t, std::shared_ptr<SweepBoundary>>& sweep_boundaries,
             bool vacuum_bcs_are_dirichlet)
{
  std::map<uint64_t, BoundaryCondition> bcs;
  for (const auto& [bid, lbs_bndry] : sweep_boundaries)
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

std::map<unsigned int, Multigroup_D_and_sigR>
PackGroupsetXS(const BlockID2XSMap& blkid_to_xs_map,
               unsigned int first_grp_index,
               unsigned int last_group_index,
               double time_absorption_scale)
{
  OpenSnInvalidArgumentIf(last_group_index < first_grp_index,
                          "last_grp_index must be >= first_grp_index");
  const unsigned int num_gs_groups = last_group_index - first_grp_index + 1;

  std::map<unsigned int, Multigroup_D_and_sigR> matid_2_mgxs_map;
  for (const auto& matid_xs_pair : blkid_to_xs_map)
  {
    const auto& mat_id = matid_xs_pair.first;
    const auto& xs = matid_xs_pair.second;

    std::vector<double> D(num_gs_groups, 0.0);
    std::vector<double> sigma_r(num_gs_groups, 0.0);

    unsigned int g = 0;
    const auto& diffusion_coeff = xs->GetDiffusionCoefficient();
    const auto& sigma_removal = xs->GetSigmaRemoval();
    const auto& inv_velocity = xs->GetInverseVelocity();
    OpenSnLogicalErrorIf(time_absorption_scale > 0.0 and inv_velocity.size() < last_group_index + 1,
                         "Time-dependent diffusion acceleration requires inverse velocities for "
                         "every material.");
    for (unsigned int gprime = first_grp_index; gprime <= last_group_index; ++gprime)
    {
      const double tau =
        time_absorption_scale > 0.0 ? inv_velocity[gprime] * time_absorption_scale : 0.0;
      D[g] = AddTimeAbsorptionToDiffusionCoefficient(diffusion_coeff[gprime], tau);
      sigma_r[g] = sigma_removal[gprime] + tau;
      ++g;
    } // for g

    matid_2_mgxs_map.insert(std::make_pair(mat_id, Multigroup_D_and_sigR{D, sigma_r}));
  }

  return matid_2_mgxs_map;
}

double
AddTimeAbsorptionToDiffusionCoefficient(const double D, const double tau)
{
  return tau > 0.0 ? D / (1.0 + 3.0 * D * tau) : D;
}

void
ApplyDSACorrectionToDelayedBoundaryFlux(DiscreteOrdinatesProblem& do_problem,
                                        const LBSGroupset& groupset,
                                        const std::vector<double>& delta_phi_local)
{
  const auto& grid = *do_problem.GetGrid();
  const auto& sdm = do_problem.GetSpatialDiscretization();
  const auto& transport_views = do_problem.GetCellTransportViews();
  const auto first_group = groupset.first_group;

  const auto correction = [&](std::uint32_t cell_local_id,
                              unsigned int face,
                              unsigned int face_node,
                              unsigned int group_idx)
  {
    const auto& cell = *grid.GetLocalCells()[cell_local_id];
    const auto node = sdm.GetCellMapping(cell).MapFaceNode(face, face_node);
    return delta_phi_local[transport_views[cell_local_id].MapDOF(node, 0, first_group + group_idx)];
  };

  for (const auto& [bid, boundary] : do_problem.GetSweepBoundaries())
    boundary->AddToNewDelayedAngularFlux(groupset.id, correction);
}

} // namespace opensn
