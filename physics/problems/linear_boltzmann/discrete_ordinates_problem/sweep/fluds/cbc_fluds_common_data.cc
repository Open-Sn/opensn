// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "physics/problems/linear_boltzmann/discrete_ordinates_problem/sweep/fluds/cbc_fluds_common_data.h"
#include "physics/problems/linear_boltzmann/discrete_ordinates_problem/sweep/spds/spds.h"
#include "framework/mesh/mesh_continuum/mesh_continuum.h"

namespace opensn
{

CBC_FLUDSCommonData::CBC_FLUDSCommonData(
  const SPDS& spds, const std::vector<CellFaceNodalMapping>& grid_nodal_mappings)
  : FLUDSCommonData(spds, grid_nodal_mappings)
{
}

} // namespace opensn
