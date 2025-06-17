// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include "physics/problems/linear_boltzmann/discrete_ordinates_problem/sweep/fluds/fluds_common_data.h"
#include <cinttypes>

namespace opensn
{

class CBC_FLUDSCommonData : public FLUDSCommonData
{
public:
  CBC_FLUDSCommonData(const SPDS& spds,
                      const std::vector<CellFaceNodalMapping>& grid_nodal_mappings);
};

} // namespace opensn
