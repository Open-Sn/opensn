// SPDX-FileCopyrightText: 2025 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/discrete_ordinates_problem.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep/fluds/aahd_fluds_common_data.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep/fluds/aahd_fluds.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep/fluds/cbcd_fluds_common_data.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep/fluds/cbcd_fluds.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep/angle_set/cbcd_angle_set.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep_chunks/cbcd_sweep_chunk.h"
#include "framework/math/spatial_discretization/spatial_discretization.h"
#include "framework/mesh/mesh_continuum/mesh_continuum.h"
#include "sweep/angle_set/cbcd_angle_set.h"

namespace opensn
{

void
DiscreteOrdinatesProblem::CreateFLUDSCommonDataForDevice()
{
  for (const auto& [quadrature, spds_list] : quadrature_spds_map_)
  {
    for (const auto& spds : spds_list)
    {
      quadrature_fluds_commondata_map_[quadrature].push_back(
        std::make_unique<AAHD_FLUDSCommonData>(*spds, grid_nodal_mappings_, *discretization_));
    }
  }
}

std::shared_ptr<FLUDS>
DiscreteOrdinatesProblem::CreateFLUDSForDevice(std::size_t num_groups,
                                               std::size_t num_angles,
                                               const FLUDSCommonData& common_data)
{
  return std::make_shared<AAHD_FLUDS>(
    num_groups, num_angles, dynamic_cast<const AAHD_FLUDSCommonData&>(common_data));
}

void
DiscreteOrdinatesProblem::CreateCBCD_FLUDSCommonData()
{
  for (const auto& [quadrature, spds_list] : quadrature_spds_map_)
  {
    for (const auto& spds : spds_list)
    {
      quadrature_fluds_commondata_map_[quadrature].push_back(
        std::make_unique<CBCD_FLUDSCommonData>(*spds, grid_nodal_mappings_, *discretization_));
    }
  }
}

std::shared_ptr<FLUDS>
DiscreteOrdinatesProblem::CreateCBCD_FLUDS(std::size_t num_groups,
                                           std::size_t num_angles,
                                           std::size_t num_local_cells,
                                           const FLUDSCommonData& common_data,
                                           const UnknownManager& psi_uk_man,
                                           const SpatialDiscretization& sdm)
{
  return std::make_shared<CBCD_FLUDS>(num_groups,
                                      num_angles,
                                      num_local_cells,
                                      dynamic_cast<const CBCD_FLUDSCommonData&>(common_data),
                                      psi_uk_man,
                                      sdm);
}

std::shared_ptr<AngleSet>
DiscreteOrdinatesProblem::CreateCBCD_AngleSet(
  size_t id,
  size_t num_groups,
  const SPDS& spds,
  std::shared_ptr<FLUDS>& fluds,
  std::vector<size_t>& angle_indices,
  std::map<uint64_t, std::shared_ptr<SweepBoundary>>& boundaries,
  const MPICommunicatorSet& in_comm_set,
  bool use_gpus)
{
  return std::make_shared<CBCD_AngleSet>(
    id, num_groups, spds, fluds, angle_indices, boundaries, in_comm_set, use_gpus);
}

std::shared_ptr<SweepChunk>
DiscreteOrdinatesProblem::CreateCBCD_SweepChunk(DiscreteOrdinatesProblem& problem,
                                                LBSGroupset& groupset)
{
  return std::make_shared<CBCD_SweepChunk>(*this, groupset);
}

} // namespace opensn
