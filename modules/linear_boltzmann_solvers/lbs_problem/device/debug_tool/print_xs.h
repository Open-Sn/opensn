// SPDX-FileCopyrightText: 2025 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include "modules/linear_boltzmann_solvers/lbs_problem/device/carrier/total_xs_carrier.h"

namespace opensn::debug
{

/// Print total cross section for each block ID on device.
void PrintTotalXSOnDevice(TotalXSCarrier& xs);

} // namespace opensn::debug
