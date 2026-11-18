// SPDX-FileCopyrightText: 2026 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep_chunks/gpu_kernel/cbcd_fused_main.h"

namespace opensn::gpu_kernel
{

template void
LaunchCBCDFusedSweepAs<1, CBCD_ALL_CELL_DOFS>(const CBCD_Launch& launch,
                                              const Arguments<SweepKind::CBC>* arguments,
                                              double* const* saved_psi,
                                              CBCD_LaunchStream& stream);

} // namespace opensn::gpu_kernel
