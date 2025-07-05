// SPDX-FileCopyrightText: 2025 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep/angle_set/angle_set.h"

namespace opensn
{

#ifndef __OPENSN_USE_CUDA__
void
AngleSet::InitializeMemoryPin()
{
}

void
AngleSet::ResetMemoryPin()
{
}
#endif // __OPENSN_USE_CUDA__

AngleSet::~AngleSet()
{
  ResetMemoryPin();
}

} // namespace opensn
