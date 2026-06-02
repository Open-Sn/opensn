// SPDX-FileCopyrightText: 2026 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep/sweep_runtime_builder.h"

#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep/fluds/aahd_fluds_common_data.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep/fluds/cbcd_fluds_common_data.h"

namespace opensn::detail
{

void
BuildAAHGPUFludsCommonData(SweepRuntime& runtime,
                           const SpatialDiscretization& discretization,
                           const std::vector<CellFaceNodalMapping>& grid_nodal_mappings)
{
  for (const auto& [quadrature, spds_list] : runtime.quadrature_spds_map)
    for (const auto& spds : spds_list)
      runtime.quadrature_fluds_commondata_map[quadrature].push_back(
        std::make_unique<AAHD_FLUDSCommonData>(*spds, grid_nodal_mappings, discretization));
}

void
BuildCBCGPUFludsCommonData(SweepRuntime& runtime,
                           const SpatialDiscretization& discretization,
                           const std::vector<CellFaceNodalMapping>& grid_nodal_mappings)
{
  for (const auto& [quadrature, spds_list] : runtime.quadrature_spds_map)
    for (const auto& spds : spds_list)
      runtime.quadrature_fluds_commondata_map[quadrature].push_back(
        std::make_unique<CBCD_FLUDSCommonData>(*spds, grid_nodal_mappings, discretization));
}

} // namespace opensn::detail
