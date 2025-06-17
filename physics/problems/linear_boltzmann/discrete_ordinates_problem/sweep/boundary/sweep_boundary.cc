// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "physics/problems/linear_boltzmann/discrete_ordinates_problem/sweep/boundary/sweep_boundary.h"
#include "framework/logging/log.h"
#include "framework/runtime.h"

namespace opensn
{

double*
SweepBoundary::PsiIncoming(uint64_t cell_local_id,
                           unsigned int face_num,
                           unsigned int fi,
                           unsigned int angle_num,
                           int group_num)
{
  throw std::runtime_error(
    "SweepBoundary: PsiIncoming call made to boundary that has no such information.");
  return nullptr;
}

double*
SweepBoundary::PsiOutgoing(uint64_t cell_local_id,
                           unsigned int face_num,
                           unsigned int fi,
                           unsigned int angle_num)
{
  throw std::runtime_error(
    "SweepBoundary: PsiOutgoing call made to boundary that has no such information.");
  return nullptr;
}

} // namespace opensn
