// SPDX-FileCopyrightText: 2025 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include "modules/linear_boltzmann_solvers/lbs_problem/device/carrier/mesh_carrier.h"

namespace opensn::device_dbg
{

/// Print mesh on device.
void PrintMeshOnDevice(MeshCarrier& mesh);

} // namespace opensn::device_dbg
