// SPDX-FileCopyrightText: 2025 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include "modules/linear_boltzmann_solvers/lbs_problem/device/carrier/quadrature_carrier.h"

namespace opensn::debug
{

/// Print quadrature on device.
void PrintQuadratureOnDevice(QuadratureCarrier& quad);

} // namespace opensn::debug
