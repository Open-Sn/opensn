// SPDX-FileCopyrightText: 2025 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/discrete_ordinates_problem.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep/fluds/aahd_fluds_common_data.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep/fluds/aahd_fluds.h"

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

} // namespace opensn
