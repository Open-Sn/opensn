// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "modules/linear_boltzmann_solvers/discrete_ordinates_solver/sweep/boundary/sweep_boundary.h"
#include "framework/logging/log.h"
#include "framework/runtime.h"

namespace opensn
{

double*
SweepBoundary::PsiIncoming(uint64_t cell_local_id,
                           unsigned int face_num,
                           unsigned int fi,
                           unsigned int angle_num,
                           int group_num,
                           size_t gs_ss_begin)
{
  log.LogAllError() << "PsiIncoming call made to boundary that has no such information.";
  Exit(EXIT_FAILURE);
  return nullptr;
}

double*
SweepBoundary::PsiOutgoing(uint64_t cell_local_id,
                           unsigned int face_num,
                           unsigned int fi,
                           unsigned int angle_num,
                           size_t gs_ss_begin)
{
  log.LogAllError() << "PsiOutgoing call made to boundary that has no such information.";
  Exit(EXIT_FAILURE);
  return nullptr;
}

} // namespace opensn
