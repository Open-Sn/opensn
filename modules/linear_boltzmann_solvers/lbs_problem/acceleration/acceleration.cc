// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "modules/linear_boltzmann_solvers/lbs_problem/acceleration/acceleration.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep/boundary/sweep_boundary.h"
#include "modules/linear_boltzmann_solvers/lbs_problem/lbs_problem.h"
#include "framework/materials/multi_group_xs/multi_group_xs.h"
#include "framework/utils/error.h"
#include "framework/runtime.h"
#include <algorithm>
#include <limits>
#include <cmath>
#include <sstream>

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
               unsigned int last_group_index)
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
    for (unsigned int gprime = first_grp_index; gprime <= last_group_index; ++gprime)
    {
      D[g] = diffusion_coeff[gprime];
      sigma_r[g] = sigma_removal[gprime];
      ++g;
    } // for g

    matid_2_mgxs_map.insert(std::make_pair(mat_id, Multigroup_D_and_sigR{D, sigma_r}));
  }

  return matid_2_mgxs_map;
}

void
CheckBlockwiseUniformDensities(const LBSProblem& lbs_problem, const std::string& acceleration_name)
{
  const auto& densities = lbs_problem.GetDensitiesLocal();
  const auto& cells = lbs_problem.GetGrid()->local_cells;
  const auto& xs_map = lbs_problem.GetBlockID2XSMap();

  std::vector<unsigned int> block_ids;
  block_ids.reserve(xs_map.size());
  std::map<unsigned int, size_t> block_id_to_index;
  size_t block_index = 0;
  for (const auto& [block_id, _] : xs_map)
  {
    block_ids.push_back(block_id);
    block_id_to_index[block_id] = block_index++;
  }

  std::vector<double> local_min(block_ids.size(), std::numeric_limits<double>::infinity());
  std::vector<double> local_max(block_ids.size(), -std::numeric_limits<double>::infinity());

  for (const auto& cell : cells)
  {
    const auto block_it = block_id_to_index.find(cell.block_id);
    if (block_it == block_id_to_index.end())
      continue;

    const auto i = block_it->second;
    const auto rho = densities[cell.local_id];
    local_min[i] = std::min(local_min[i], rho);
    local_max[i] = std::max(local_max[i], rho);
  }

  std::vector<double> global_min(block_ids.size(), 0.0);
  std::vector<double> global_max(block_ids.size(), 0.0);
  OpenSnInvalidArgumentIf(block_ids.size() > static_cast<size_t>(std::numeric_limits<int>::max()),
                          "Block count exceeds supported MPI all_reduce count range.");
  const auto num_blocks = static_cast<int>(block_ids.size());
  mpi_comm.all_reduce(local_min.data(), num_blocks, global_min.data(), mpi::op::min<double>());
  mpi_comm.all_reduce(local_max.data(), num_blocks, global_max.data(), mpi::op::max<double>());

  for (size_t i = 0; i < block_ids.size(); ++i)
  {
    if ((not std::isfinite(global_min[i])) or (not std::isfinite(global_max[i])))
      continue;

    const auto tol = 1.0e-12 * std::max({1.0, std::fabs(global_min[i]), std::fabs(global_max[i])});
    if (std::fabs(global_max[i] - global_min[i]) > tol)
    {
      std::stringstream msg;
      msg << acceleration_name << " requires uniform density per block ID. Block " << block_ids[i]
          << " has non-uniform density range [" << global_min[i] << ", " << global_max[i] << "].";
      OpenSnInvalidArgument(msg.str());
    }
  }
}

} // namespace opensn
