// SPDX-FileCopyrightText: 2026 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep/boundary/arbitrary_boundary.h"
#include "modules/linear_boltzmann_solvers/lbs_problem/groupset/lbs_groupset.h"
#include <utility>

namespace opensn
{

ArbitraryBoundary::ArbitraryBoundary(BoundaryBank& bank,
                                     const std::vector<LBSGroupset>& groupsets,
                                     std::shared_ptr<AngularFluxFunction> angular_flux_function)
  : SweepBoundary(bank, LBSBoundaryType::ARBITRARY),
    angular_flux_function_(std::move(angular_flux_function))
{
  extra_data_.reserve(groupsets.size());
  extra_data_.resize(groupsets.size());

  for (const auto& groupset : groupsets)
  {
    auto quadrature = groupset.quadrature;
    auto num_angles = quadrature->omegas.size();

    bank_.ExtendBoundaryFlux(groupset.id, num_angles * groupset.GetNumGroups());
    auto& map_dirnum = extra_data_[groupset.id].map_dirnum;
    map_dirnum.reserve(num_angles);
    map_dirnum.assign(num_angles, std::numeric_limits<std::uint64_t>::max());

    auto& counter = bank_[groupset.id].counter;
    offset_[groupset.id] = counter;
    counter += num_angles;
  }
}

void
ArbitraryBoundary::InitializeAngleDependent(const std::vector<LBSGroupset>& groupsets)
{
  for (const auto& groupset : groupsets)
  {
    const auto& quadrature = groupset.quadrature;
    auto num_angles = quadrature->omegas.size();
    auto& map_dirnum = extra_data_[groupset.id].map_dirnum;
    double* boundary_flux = GetBoundaryFlux(groupset.id);

    auto angle_agg = groupset.angle_agg;
    for (std::uint64_t i = 0; const auto& angleset : *angle_agg)
    {
      const auto& angle_indices = angleset->GetAngleIndices();
      for (const auto& angle_indice : angle_indices)
        map_dirnum[angle_indice] = i++;
    }

    for (std::uint64_t dir_num = 0; dir_num < num_angles; ++dir_num)
    {
      auto internal_angle_index = map_dirnum[dir_num];
      double* angular_boundary_flux =
        boundary_flux + internal_angle_index * groupset.GetNumGroups();
      for (unsigned int g = 0; g < groupset.GetNumGroups(); ++g)
      {
        unsigned int group = groupset.first_group + g;
        angular_boundary_flux[g] =
          (*angular_flux_function_)(static_cast<int>(group), static_cast<int>(dir_num));
      }
    }
  }
}

double*
ArbitraryBoundary::PsiIncoming(std::uint32_t cell_local_id,
                               unsigned int face_num,
                               unsigned int fi,
                               unsigned int angle_num,
                               int groupset_id,
                               unsigned int group_idx)
{
  auto& extra_data = extra_data_[groupset_id];
  auto internal_angle_index = extra_data.map_dirnum[angle_num];
  return GetBoundaryFlux(groupset_id, internal_angle_index) + group_idx;
}

std::uint64_t
ArbitraryBoundary::GetOffsetToAngleset(const FaceNode& face_node,
                                       AngleSet& anglset,
                                       bool is_outgoing)
{
  if (is_outgoing)
    throw std::logic_error("ArbitraryBoundary does not support outgoing flux.");
  int groupset_id = anglset.GetGroupsetID();
  auto& extra_data = extra_data_[groupset_id];
  auto internal_angle_index = extra_data.map_dirnum[anglset.GetAngleIndices()[0]];
  return (offset_[groupset_id] + internal_angle_index) * bank_[groupset_id].groupset_size;
}

} // namespace opensn
